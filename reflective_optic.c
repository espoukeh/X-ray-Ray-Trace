/*
* reflective_optic
NOTE: NO phase change upon reflection included yet!!!!

*/


#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "ray.h"

#define PI 3.14159265358979323846264338327

static double tol = 1.0E-12;
static double limit = 10.0;
static int MAX_ITER = 100;

static int compare_fabs (const void *a, const void *b)
{
	//comparison function for qsort
	//note use of absolute value. 
	return ( fabs(*(double*)a) - fabs(*(double*)b) );	
}

static double function_eval(double *a, double *b, double *c, double *d, int N, Vector v)
{
/*
*  Evaluates function and returns result
*  Attempts are made to reduce roundoff errors by sorting the values of the partial results
*  Then adding them from the absolute value of the smallest to the largest.
*/
	double *partial_result;
	double result = 0.0;
	double compensation = 0.0;
	int i;
	partial_result  = (double *) calloc(N, sizeof(double));
	for (i = 0; i < N; i++)	
		partial_result[i] = (double) (a[i] * (double)pow(v.x,b[i]) * (double)pow(v.y, c[i]) * (double)pow(v.z, d[i]));

	//sorts the partial results into smallest to largest (using absolute values)
	//this is so that the sum is independent of the term order of the equation
	if (N > 2)
		qsort(partial_result, N, sizeof(double), compare_fabs); 

	// use of Kahan summation algorithm to reduce roundoff error	
	for (i = 0; i < N; i++) {
		double y = partial_result[i] - compensation;
		double t = result + y;
		compensation = (t - result) - y;
		result = t;
	}

	free(partial_result);
	return result;

}

/*
static double derivative_eval(double *a, double *b, double *c, double *d, int N, Vector v, Vector u)
{

//  Evaluates derivative of function and returns result
//  Attempts are made to reduce roundoff errors by sorting the values of the partial results
//  Then adding them from the absolute value of the smallest to the largest.

	double *partial_result;
	double result = 0.0;
	int i;
	partial_result  = (double *) calloc(N, sizeof(double));
	for (i = 0; i < N; i++)	
		partial_result[i] = a[i]*((b[i]*u.x*pow(v.x,b[i]-1.0)*pow(v.y,c[i])*pow(v.z,d[i]))+
		                          (c[i]*u.y*pow(v.x,b[i])*pow(v.y,c[i]-1.0)*pow(v.z,d[i]))+
		                          (d[i]*u.z*pow(v.x,b[i])*pow(v.y,c[i])*pow(v.z,d[i]-1.0)));

	//sorts the partial results into smallest to largest (using absolute values)
	if (N > 2)
		qsort(partial_result, N, sizeof(double), compare_abs); 
	for (i = 0; i < N; i++)
		result += partial_result[i];
	free(partial_result);
	return result;
}
*/

static double find_surface(Ray R, double *a, double *b, double *c, double *d, int N)
{
/* 
* Returns the propagation distance,L, to the surface of the optic.
* If there is an error, it returns NAN.
* Uses bisection method. If (i,j,k) is the current ray position and (u,v,w)
* is the direction of the ray. All points, (x,y,z), allong the ray can be expressed as:
* (x,y,z) = (i+(L*u), j+(L*v), k+(L*w)). Where L is a real value.
* Therefor function f(x,y,z) = 0 is a 1D function of L where the root is the
* intersection with the optical surface.
* 
* Note: Newton method was not used as it occassionally will propagated backwards 
* (all real optics have 2 solutions).
*/
	double Lmin, Lmax, Lmid;
	double Vmin, Vmax, Vmid;
	Vector V;
	int i;
	
	Lmin = 0.0;
	V = make_vector(R.pos.x + (Lmin * R.v.x),
                        R.pos.y + (Lmin * R.v.y),
                        R.pos.z + (Lmin * R.v.z));
	Vmin = 	function_eval(a,b,c,d,N,V);

	Lmax = limit;
	V = make_vector(R.pos.x + (Lmax * R.v.x),
                        R.pos.y + (Lmax * R.v.y),
                        R.pos.z + (Lmax * R.v.z));
	Vmax = function_eval(a,b,c,d,N,V);
	if(Vmin * Vmax > 0.0) return NAN; //function limits are wrong!

	for (i = 0; i < MAX_ITER; i++)
	{		
		Lmid = (Lmax + Lmin) / 2.0;
		Vector V = make_vector(R.pos.x + (Lmid * R.v.x),
		                       R.pos.y + (Lmid * R.v.y),
		                       R.pos.z + (Lmid * R.v.z));

		Vmid = function_eval(a,b,c,d,N,V);

		if (fabs(Vmid) < tol) return Lmid;
		if (Vmid * Vmax > 0.0) {
			Lmax = Lmid;
			Vmax = Vmid;
		} else {
			Lmin = Lmid;
			Vmin = Vmid;
		}
	}
	return NAN; //maximum iteration reached!
}

static int read_function(char *function, double *a, double *b, double *c, double *d, int N)
{
/* Reads function from command line assuming regular format */
	int rest;
	int i,j;
	fprintf(stderr,"Reading function:\n");
	for (i = 0; i < N; i++) {
		j = sscanf(function, "%lf*x^%lf*y^%lf*z^%lf %n",&(a[i]),&(b[i]),&(c[i]),&(d[i]),&rest);
		if (i != N-1) {
			fprintf(stderr,"%f*x^%f*y^%f*z^%f + ",a[i],b[i],c[i],d[i]);
		} else {
			fprintf(stderr,"%f*x^%f*y^%f*z^%f\n",a[i],b[i],c[i],d[i]);
		}
		if (j != 4) {
			fprintf(stderr,"Function formatting error!\n"
			               "could not indentify function: %s\n", function);
			return -1;
		}
		function += rest;
	}
	return 0;
}

static int make_partial_derivative(double *a, double *b, double *dxa, double *dxb, int N)
{
/* 
*  Makes partial derivatives
*  If f(x) = a[0]*x^b[0] + a[1]*x^b[1] + ... + a[n-1]*x^b[n-1]
*  Then f'(x) = dxa[0]*x^dxb[] + a[1]*x^dxb[1] + ... + dxa[n-1]*x^dxb[n-1]
*  Where dxa[k] = a[k]*b[k] and dxb = b[k]-1.0
*/

	int i;
	for (i = 0; i < N; i++) {
//if f(x) = a*x^0 then f'(x) = 0, I get f'(x) = 0*x^-1; if x=0 this gives NaN!
//Best be explicit as this is a common case!
		if (b[i] != 0.0) {
			dxa[i] = a[i] * b[i];
			dxb[i] = b[i] - 1.0;
		} else {
			dxa[i] = 0.0;
			dxb[i] = 0.0;
		}
	}
	return 0;
}

static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Calculates reflected rays from arbitrary optical surfaces. The optical surface\n"
"is approximated by a polynomial function, <fxn>, of x,y,z such that: \n"
"\n"
"f(x,y,z)= 0 = a_0*x^b_0*y^c_0*z^d_0 +a_1*x^b_1*y^c_1*z^d_1 +\n"
"               ... +a_N*x^b_N*y^c_N*z^d_N\n"
"\n"
"Where N is the order of the function to be evaluate. The surface intersect\n"
"of the ray is founds using a simple bisecton method. While the normal is \n"
"determined via the gradiant of the function. Translation: Please choose your\n"
"function and strating rays such that the results exist and are physical.\n"
"Note: All values are metric.\n"
"\n"
"  -h, --help               Display this help message.\n"
"  -i, --input=<file>       Input filename. Default: stdin.\n"
"  -o, --output=<file>      Output filename. Default: stdout.\n"
"  -n, --terms=<int>        Number of terms, N, in function <fxn>\n"
"  -f, --fxn='<fxn>'        The function the describes the optics surface.\n"
"  -l, --limit<num>         The limit for propagation for bisection: the program\n"
"                             propagate the Ray from [0, l) looking for \n"
"                             intersection with fxn (default = 10.0 m)\n"
"  -t, --tolerance=<num>    Tolerance value of the optic surface. It should be a\n"
"                             fraction of the wavelength (default = 1.0E-12 m)\n"
"  -m, --iterations=<int>   Max number of iterations for bisection method \n"
"                             (default = 100)\n"
"  -r, --reflectivity=<num> Reflectivty of optics surface (defalut = 1.0)\n" 
"  -z, --zero               Enables propagation of 0 intensity rays\n"
"                             By default not propaged.\n" 
"\n");
}

int main(int argc, char *argv[])
{
	/* initial variable */
	FILE *in = stdin;
	FILE *out = stdout;
	char *inFileName = NULL;
	char *outFileName = NULL;
	char *function = NULL;
	int g;
	int N = 0;
	int useZero = 0;
	double ref = 1.0;
	double *a, *b, *c, *d; //array of values for function evaluation
	double *dxa, *dxb, *dya, *dyc, *dza, *dzd; //Partial derivative variables
	
	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"terms",             1, NULL,               'n'},
		{"fxn",               1, NULL,               'f'},
		{"limit",             1, NULL,               'l'},
		{"iterations",        1, NULL,               'm'},
		{"tolerance",         1, NULL,               't'},
		{"reflectivity",      1, NULL,               'r'},
		{"zero",              0, NULL,               'z'},
		{0, 0, NULL, 0}
	};


	 /* Short options */
        while ((g = getopt_long(argc, argv, "hi:o:n:f:l:m:t:r:z",
                longopts, NULL)) != -1)
        {

                switch (g) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'i' :
			inFileName = strdup(optarg);
			break;

			case 'o' :
			outFileName = strdup(optarg);
			break;

			case 'n' :
			N = atoi(optarg);
			if (N < 0) {
				fprintf(stderr, "Number of terms in function "
				   "needs to be at least 1.\n");
				show_help(argv[0]);
				return 1;
			}
			break;

			case 'm' :
			MAX_ITER = atoi(optarg);
			if (MAX_ITER < 1)
				fprintf(stderr, "Number of iterations"
				   "needs to be at least 1.\n");
				show_help(argv[0]);
				return 1;
			break;

			case 'f' :
			function = strdup(optarg);
			break;

			case 'l' :
			limit = atof(optarg);
			break;

			case 't' :
			tol = atof(optarg);
			break;

			case 'r' :
			ref = atof(optarg);
			if (ref > 1.0 || ref <0.0) {
				fprintf(stderr,"Invalid reflectivity: '%f'\n", ref);
				return 1;
			}
			break;

			case 'z' :
			useZero = 1;
			break;

			default :
                        fprintf(stderr,"Unexpected arguement \n");
                        show_help(argv[0]);
                        return 1;
		}
	}

	/* read function into memory & compute grad(f(x,y,z))*/
	if (function == NULL) {
		fprintf(stderr, "No optical funtion has been inputed!\n");
		show_help(argv[0]);
		return 1;
	}

	//N += 1; //add term 0
	a   = (double *) calloc(N, sizeof(double));
	b   = (double *) calloc(N, sizeof(double));
	c   = (double *) calloc(N, sizeof(double));
	d   = (double *) calloc(N, sizeof(double));
	dxa = (double *) calloc(N, sizeof(double));
	dxb = (double *) calloc(N, sizeof(double));
	dya = (double *) calloc(N, sizeof(double));
	dyc = (double *) calloc(N, sizeof(double));
	dza = (double *) calloc(N, sizeof(double));
	dzd = (double *) calloc(N, sizeof(double));
	if (a == NULL || b == NULL || c == NULL || d == NULL || dxa == NULL  ||
	    dxb == NULL || dya == NULL || dyc == NULL || dza == NULL || dzd == NULL) {
		fprintf(stderr, "Out of memory error! \n");
		return 1;
	}
	g = read_function(function, a, b, c, d, N);
	if (g != 0) return 1; //Error reading function, perhaps I should free memory?
	make_partial_derivative(a, b, dxa, dxb, N); //partial of x
	make_partial_derivative(a, c, dya, dyc, N); //partial of y
	make_partial_derivative(a, d, dza, dzd, N); //partial of z

	/* open files if necessary */
	if ( inFileName != NULL ) {
		in = fopen(inFileName, "rb");
		if ( in == NULL ) {
			fprintf(stderr,"Failed to open input file: %s!\n", inFileName);
			return 1;
		}
		free(inFileName);
	}

	if ( outFileName != NULL ) {
		out = fopen(outFileName, "wb");
		if ( out == NULL ) {
			fprintf(stderr,"Failed to open output file: %s!\n", outFileName);
			return 1;
		}
		free(outFileName);
	}

	/* loop over vectors */
	while (!feof(in)) {
		double prop, j, k, theta, VdotN;
		Ray myRay;
		Vector Normal, Out;
		g = read_ray(in, &myRay);
		if (g != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}

		/* check if using 0 or skip */
		if (useZero == 0 && myRay.i == 0.0) {
			write_ray(out, myRay);
			continue;
		}

		myRay.v = unit(myRay.v); //check that direction is a unit vector
		
		//propagate to surface
		prop = find_surface(myRay, a, b, c, d, N);

		// error check
		if(!isfinite(prop)) {
			fprintf(stderr,"Warning: can't find mirror position "
			  "for Ray %lu!\n Skipping & setting intensity to 0\n",
		           myRay.ray);
			myRay.i = 0.0;
			write_ray(out, myRay);
			continue;
		}

		//update propagation
		myRay.pos = make_vector(myRay.pos.x + (prop * myRay.v.x),
		                        myRay.pos.y + (prop * myRay.v.y),
		                        myRay.pos.z + (prop * myRay.v.z));
		myRay.l += prop;
		myRay.i = myRay.i * ref;

		//make normal vector (gradiant)
		Normal.x = function_eval(dxa, dxb, c, d, N, myRay.pos);
		Normal.y = function_eval(dya, b, dyc, d, N, myRay.pos);
		Normal.z = function_eval(dza, b, c, dzd, N, myRay.pos);
		Normal = unit(Normal);

		myRay.v = mult(myRay.v, -1.0);  //Note math easier if I flip my incident ray
		if(acos(dot(Normal,myRay.v)) > PI/2.0) 
			Normal = mult(Normal, -1.0); //Normal facing wrong way

		theta = acos(dot(Normal, myRay.v));
		VdotN = dot(myRay.v, Normal); //note both unit vectors
		if (VdotN == 1.0) {
			//vectors are colinear! Normal reflection!
			k = 1;
			j = 0;
		} else {
			k = ((cos(2.0*theta)*VdotN)-cos(theta))/(pow(VdotN,2.0) - 1.0);
			j = cos(2.0*theta) - (k*VdotN);
		}
		Out = add(mult(myRay.v,j),mult(Normal,k));
		Out = unit(Out); // just checking, but needless otherwise math is wrong
		myRay.v = Out;
		
		//print ray
		write_ray(out, myRay);
	}

	fclose(in);
	fclose(out);		
	free(a); free(b); free(c); free(d); 
	free(dxa); free(dxb); free(dya);
	free(dyc); free(dza); free(dzd);
	return 0;
}

