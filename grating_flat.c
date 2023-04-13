/*
* grating_flat
NOTE: Check phase change upon reflection!!!!
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "ray.h"
#define PI 3.14159265358979323846264338327

static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Creates a rectangular refletive grating with a variable d spacing.\n"
"Diffracts all rays in grating aperture with order m, and it sets all\n"
"rays outside aperture to 0 intensity. The grating is defined by a\n"
"rectangle with with vertices P0, P1, P2, and (P0 + P1 + P2). With\n"
"P1-P0 the x direction and P2-P0 the y direction.\n"
"The grating d spacing and direction are defined by a polynomial \n"
"function, <fxn>, of x,y such that: \n"
"\n"
"d(x,y) = a_0*x^b_0*y^c_0 + a_1*x^b_1*y^c_1 + ...\n"
"   ... + a_N*x^b_N*y^c_N\n"
"where 0<x<1 and 0<y<1 and with the origin corresponding to P0 \n"
"and (1,1) corresponding to the x,y positions of (P0 + P1 + P2).\n"
"Grating direction is in the GRAD(d(x,y)) direction. Program does\n"
"NOT work with constant D-spacing gratings."
"Note: all values are metric.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"      --P0='[x,y,z]'      First vertex (default = [0,0,0]).\n"
"      --P1='[x,y,z]'      Second vertex (default = [1,0,0]).\n"
"      --p2='[x,y,z]'      Third vertex (default = [0,1,0]).\n"
"  -n, --terms=<int>       Number of terms, N, in d-spacing <fxn>\n"
"  -d, --dspacing=<fxn>    D-spacing [m].\n"
"  -m, --order=<int>       Grating order (default = 1)\n"
"  -z, --zero              Enables propagation of 0 intensity rays\n"
"                            By default not propaged.\n" 
"\n");
}

static int read_function(char *function, double *a, double *b, double *c, int N)
{
/* Reads function from command line assuming regular format */
	int rest;
	int i,j;
	fprintf(stderr,"Reading function:\n");
	for (i = 0; i < N; i++) {
		j = sscanf(function, "%lf*x^%lf*y^%lf %n",&(a[i]),&(b[i]),&(c[i]),&rest);
		if (i != N-1) {
			fprintf(stderr,"%f*x^%f*y^%f + ",a[i],b[i],c[i]);
		} else {
			fprintf(stderr,"%f*x^%f*y^%f\n",a[i],b[i],c[i]);
		}
		if (j != 3) {
			fprintf(stderr,"Function formatting error!\n"
			               "could not indentify function: %s\n", function);
			return -1;
		}
		function += rest;
	}

	//simple sanity check thata a grating direction exists.
	j = 0;
	for (i = 0; i < N; i++) {
		if (b[i] != 0.0 || c[i] != 0.0) {
			j = 1;
			break;
		}
	}
	if (j == 0) {
		fprintf(stderr,"D-spacing function is a constant.\n");
		return -1;
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

static int compare_fabs (const void *a, const void *b)
{
	//comparison function for qsort
	//note use of absolute value. 
	return ( fabs(*(double*)a) - fabs(*(double*)b) );	
}

static double function_eval(double *a, double *b, double *c, int N, double x, double y)
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
		partial_result[i] = (double) (a[i] * (double)pow(x,b[i]) * (double)pow(y, c[i]));

	//sorts the partial results into smallest to largest (using absolute values)
	//this is so that the sum is independent of the term order of the equation
	if (N > 2)
		qsort(partial_result, N, sizeof(double), compare_fabs); 

	// use of Kahan summation algorithm to reduce roundoff error	
	for (i = 0; i < N; i++) {
		double s = partial_result[i] - compensation;
		double t = result + s;
		compensation = (t - result) - s;
		result = t;
	}

	free(partial_result);
	return result;
}

int main(int argc, char *argv[])
{
	/* initial variable */
	//file stuff
	FILE *in = stdin;
	FILE *out = stdout;
	char *inFileName = NULL;
	char *outFileName = NULL;

	// input options stuff
	char *function = NULL;
	int g;
	int N = 0;
	int useZero = 0;
	int order = 1;
	double *a, *b, *c; //array of values for function evaluation
	double *dxa, *dxb, *dya, *dyc; //Partial derivative variables
	Vector P0 = make_vector(0.0, 0.0, 0.0);
	Vector P1 = make_vector(1.0, 0.0, 0.0);
	Vector P2 = make_vector(0.0, 1.0, 0.0);

	// math stuff
	Vector B1, B2, B1unitV, B2unitV, Normal; //basis vectors & normal
	double B1mag, B2mag;

	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"P0",                1, NULL,               'a'},
		{"P1",                1, NULL,               'b'},
		{"P2",                1, NULL,               'c'},
		{"dspacing",          1, NULL,               'd'},
		{"terms",             1, NULL,               'n'},
		{"order",             1, NULL,               'm'},
		{"zero",              0, NULL,               'z'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((g = getopt_long(argc, argv, "hi:o:a:b:c:d:n:m:z",
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

			case 'a' :
			sscanf(optarg,"[%le,%le,%le]", &P0.x, &P0.y, &P0.z);
			break;

			case 'b' :
			sscanf(optarg,"[%le,%le,%le]", &P1.x, &P1.y, &P1.z);
			break;

			case 'c' :
			sscanf(optarg,"[%le,%le,%le]", &P2.x, &P2.y, &P2.z);
			break;

			case 'd' :
			function = strdup(optarg);
			break;

			case 'm' :
			order = atoi(optarg);
			break;

			case 'n' :
			N = atoi(optarg);
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
	
	/* read function into memory & compute grad(f(x,y))*/
	if (function == NULL) {
		fprintf(stderr, "No optical funtion has been inputed!\n");
		show_help(argv[0]);
		return 1;
	}

	if (N < 0) {
		fprintf(stderr, "Number of terms in function "
		   "needs to be at least 1.\n");
		show_help(argv[0]);
		return 1;
	}
	//N += 1; //add term 0

	a   = (double *) calloc(N, sizeof(double));
	b   = (double *) calloc(N, sizeof(double));
	c   = (double *) calloc(N, sizeof(double));
	dxa = (double *) calloc(N, sizeof(double));
	dxb = (double *) calloc(N, sizeof(double));
	dya = (double *) calloc(N, sizeof(double));
	dyc = (double *) calloc(N, sizeof(double));
	if (a == NULL || b == NULL || c == NULL || dxa == NULL  ||
	    dxb == NULL || dya == NULL || dyc == NULL) {
		fprintf(stderr, "Out of memory error! \n");
		return 1;
	}

	g = read_function(function, a, b, c, N);
	if (g != 0) return 1; //Error reading function, perhaps I should free memory?
	make_partial_derivative(a, b, dxa, dxb, N); //partial of x
	make_partial_derivative(a, c, dya, dyc, N); //partial of y

	//calculate cross product sanity check
	Normal = cross(sub(P1,P0),sub(P2,P0));

	if (magnitude(Normal) == 0.0) {  //see if points are colinear
		fprintf(stderr,"Error points are colinear!\n");
		return 1;
	}
	Normal = unit(Normal);

	if (dot(sub(P1,P0),sub(P2,P0)) != 0.0) { //points are not perpendicular
		fprintf(stderr,"Error points are not perpendicular!\n");
		// NOTE: suggest direction to set P2?
		return 1;
	}

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

	/* calculate basis vectors in grating plane*/
	B1 = sub(P1,P0);
	B1unitV = unit(B1);
	B1mag = magnitude(B1);
	B2 = sub(P2,P0);
	B2unitV = unit(B2);
	B2mag = magnitude(B2);
	
	/* loop over vectors */
	while (!feof(in)) {				

		double distance, x, y;
		Ray myRay;
		Vector point;

		g = read_ray(in, &myRay);
		if (g != 0) {
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}
		
		/* check if using 0 intensity or skip */
		if (useZero == 0 && myRay.i == 0.0) {
			write_ray(out, myRay);
			continue;
		}

		/* check ray is not parallel to mirror */
		if (dot(myRay.v,Normal) == 0.0) {
			fprintf(stderr,"Warning: Ray %lu propagates parallel to"
			 " surface!\n Skipping ray, and copying it to output.\n", myRay.ray);
			write_ray(out, myRay);
			continue;
		} 

		myRay.v = unit(myRay.v); //make the direction a unit vector
		//caluclate distance from the location point in myRay to 
		//grating plane along propagation direction
		distance = dot(sub(P0,myRay.pos),Normal)/dot(myRay.v,Normal);
		if (distance < 0.0) {
			fprintf(stderr,"Warning Ray %lu is furthern then the grating's surface"
			 " back propagating to grating.\n", myRay.ray);
		}
		myRay.pos = add(myRay.pos,mult(myRay.v,distance));
		myRay.l = myRay.l + distance;

		//calculate where myRay.pos is with respect to P0 & basis vectors
		//B1=(P1-P0), B2=(P2-P0) & dot(B1,B2)=0.
		point = sub(myRay.pos, P0);
		x = dot(point,B1unitV)/B1mag; //allowed as dot(B1,B2)=0
		y = dot(point,B2unitV)/B2mag;	
		if ( x >= 0.0 && x <= 1.0 && y >= 0.0 && y <= 1.0 ) {
			//in grating
			double i,j,k, theta, alpha, beta, dSpacing, dx, dy;
			Vector Out, G, D;

			//calculate d-spacing		
			dSpacing = function_eval(a, b, c, N, x, y);

			//calculate direction of grating
			//Note see if using U.W.Ludwig JOSA 63, 1105 (1973) would be better

			dx = function_eval(dxa, dxb, c, N, x, y);
			dy = function_eval(dya, b, dyc, N, x, y);

			//direction check
			if (dx == 0.0 && dy == 0.0) {
			fprintf(stderr,"Warning can't find grating direction for Ray: %lu \n"
			 "Skipping ray & setting intensity to 0, but copying it to output.\n", myRay.ray);
				myRay.i = 0.0;				
				write_ray(out, myRay);
				continue;
			}

			G = unit(add(P0, add( mult(B1, dx), mult(B2, dy)))); //unit vector in grating plane, perp to constant grating rulling
			D = unit(cross(Normal, G)); //unit vector parallel to grating & in plane

			myRay.v = mult(myRay.v, -1.0);  //Note math easier if I flip my incident ray
			if(acos(dot(Normal,myRay.v)) > PI/2.0) 
				Normal = mult(Normal, -1.0); //Normal facing wrong way

			//Calculate angles			
			theta = atan(dot(G, myRay.v) / dot(Normal, myRay.v)); //incidence angle perp to grating
			alpha = atan(dot(D, myRay.v) / dot(Normal, myRay.v)); //incidence angle parallel to grating
			beta = asin((order*myRay.w/dSpacing) - sin(alpha)); //diffraction angle

			i = cos((PI/2.0)-beta);
			j = cos((PI/2.0)-theta);
			k = sqrt(1.0 - (i*i) - (j*j)); 
			//Note there are two solutions (-k works too)! 
			//I am taking the positive "reflective" one as Normal is defined outside surface is positive

			// i*B1unitV + j*B2unitV + k*Normal = Out_vector
			Out = add(mult(B1unitV,i),mult(B2unitV,j));
			Out = add(Out,mult(Normal, k));
			Out = unit(Out); // just checking, but needless if math above is correct
			myRay.v = Out;
//fprintf(stderr,"alpha = %f, beta = %f, theta = %f\n", alpha*180/PI, beta*180/PI, theta*180/PI);
//fprintf(stderr,"i = %f, j = %f, k = %f\n",i,j,k);
//fprintf(stderr,"Out = (%f, %f, %f)\n",Out.x, Out.y, Out.z);
		} else {
			myRay.i = 0.0; // ray is outside mirror aperture
		}

		//print ray
		write_ray(out, myRay);
	}

	free(a); free(b); free(c); 
	free(dxa); free(dxb);
	free(dya); free(dyc);
	fclose(in);
	fclose(out);
	return 0;
}
