/*
* reflective_shape
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

enum {
	NONE,
	LINEAR,
	CUBIC
};

static int compare_fabs (const void *a, const void *b)
{
	//comparison function for qsort
	//note use of absolute value. 
	return ( fabs(*(double*)a) - fabs(*(double*)b) );	
}

static int read_mirror(double **profile, char *name, int rows, int columns)
{
/* Reads mirror profile 
   NOTE: little error checking is done on the profile should fix this */
	FILE *mirror = NULL;
	double value;
	int i,j;
	mirror = fopen(name, "r");
	if ( mirror == NULL ) {
		fprintf(stderr,"Failed to open input file: %s!\n", name);
		return 1;
	}
	
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++) {
			fscanf(mirror, " %lf , ", &value);
			profile[i][j] = value;
		}
	}
	fclose(mirror);
	free(name);
	return 0;
}

static double max_mirror(double **m, int rows, int columns)
{
	double max;
	int i, j;
	max = m[0][0];
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++) {
			 if (m[i][j] > max) {
				max = m[i][j];
			}
		}
	}
	return max;
}

static double min_mirror(double **m, int rows, int columns)
{
	double min;
	int i, j;
	min = m[0][0];
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++) {
			 if (m[i][j] < min) {
				min = m[i][j];
			}
		}
	}
	return min;
}


static void calc_b00(double **b00, double **mZ, int rows, int columns)
{
	int i,j;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++) {
			b00[i][j]= mZ[i][j];
		}
	}
}
static void calc_b01(double **b01, double **mZ, int rows, int columns)
{
//note using calloc and also don't care about boundries that won't be accessed
	int i,j;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns-1; j++) {
			b01[i][j]= mZ[i][j+1] - mZ[i][j];
		}
	}
}
static void calc_b10(double **b10, double **mZ, int rows, int columns)
{
	int i,j;
	for (i = 0; i < rows-1; i++) {
		for (j = 0; j < columns; j++) {
			b01[i][j]= mZ[i+1][j] - mZ[i][j];
		}
	}
}
static void calc_b11(double **b11, double **mZ, int rows, int columns)
{
	int i,j;
	for (i = 0; i < rows-1; i++) {
		for (j = 0; j < columns-1; j++) {
			b11[i][j]= mZ[i][j] - mZ[i+1][j] - mZ[i][j+1] + mZ[i+1][j+1];
		}
	}
}


static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Calculates reflected rays from from a surface profile from a file. The optical surface\n"
"is approximated by a n row by m columns CSV file giving the optics heigth as a\n"
"functon of position. The values in between the points are interpolated. The optic is\n"
"loacated at P0, P1, P2, and ( P0 + P1 + P2) and consists of n pixels along the \n"
"(P1-P0) direction and m pixels along the (P2-P0) direction. Of importants is\n"
"that the direction (P1-P0) is perpendicular to (P2-P0), and a right hand rule indicates\n"
"the + direction for the surface height. Also the point P0 is the first value in the\n"
"CSV file, P1 is the first row last column, and P2 is the first column last row.\n"
"The exact optic might deviate from these positions if the heigth at these points is not 0.\n"
"Note: All values are metric.\n"
"\n"
"  -h, --help               Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"      --mirror=<file>     CSV file for optic.\n"
"      --interp=<type>     Interpolation type to use between mirror points\n"
"                          Choose from:\n"
"                            none        : Nearest-neighbor interpolation\n"
"                            linear      : bilinear interpolation\n"
"                            cubic       : bicubic interpolation (default)\n"
"      --P0='[x,y,z]'      First vertex (default = [0,0,0].\n"
"      --P1='[x,y,z]'      Second vertex (default = [1,0,0].\n"
"      --p2='[x,y,z]'      Third vertex (default = [0,1,0].\n"
"  -n, --rows=<int>        Number of rows in CSV \n"
"                           (along direction between P0 and P1)\n"
"  -m, --columns=<int>     Number of columns in CSV  \n"
"                           (along direction between P0 and P2)\n"
"  -t, --tolerance=<num>    Tolerance value of the optic surface. It should be a\n"
"                             fraction of the wavelength (default = 1.0E-12 m)\n"
"      --iterations=<int>   Max number of iterations for bisection method \n"
"                             (default = 100)\n"
"  -r, --reflectivity=<num> Reflectivty of optics surface (defalut = 1.0)\n" 
"  -z, --zero               Enables propagation of 0 intensity rays\n"
"                             By default not propaged.\n");
}

int main(int argc, char *argv[])
{
	/* initial variable */
	FILE *in = stdin;
	FILE *out = stdout;
	char *inFileName = NULL;
	char *outFileName = NULL;
	char *mirrorFileName = NULL;
	double **mZ; //Mirror's positions
//	double **a00, **a01, **a02, **a03, **a10, **a11, **a12, **a12; //Use for bicubic
//	double **a20, **a21, **a22, **a23, **a30, **a31, **a32, **a33;
	double **b00, **b01, **b10, **b11; //use for bilinear
	int g;
	int rows = 0;
	int columns = 0;
	Vector P0 = make_vector(0.0, 0.0, 0.0);
	Vector P1 = make_vector(1.0, 0.0, 0.0);
	Vector P2 = make_vector(0.0, 1.0, 0.0);
	Vector B0, B1, B2, B0Unit, B1Unit, B2Unit, P3, P4;
	int useZero = 0;
	double ref = 1.0;
	double mirrorMax, mirrorMin;
	
	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"mirror",            1, NULL,               'p'},
		{"rows",              1, NULL,               'n'},
		{"columns",           1, NULL,               'm'},
		{"limit",             1, NULL,               'l'},
		{"iterations",        1, NULL,               's'},
		{"tolerance",         1, NULL,               't'},
		{"reflectivity",      1, NULL,               'r'},
		{"zero",              0, NULL,               'z'},
		{0, 0, NULL, 0}
	};


	 /* Short options */
        while ((g = getopt_long(argc, argv, "hi:o:p:n:m:l:s:t:r:z",
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

			case 'p' :
			mirrorFileName = strdup(optarg);
			break;

			case 'n' :
			rows = atoi(optarg);
			break;

			case 'm' :
			columns = atoi(optarg);
			break;

			case 's' :
			MAX_ITER = atoi(optarg);
			if (MAX_ITER < 1)
				fprintf(stderr, "Number of iterations"
				   "needs to be at least 1.\n");
				show_help(argv[0]);
				return 1;
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

			case 'z' :
			useZero = 1;
			break;

			default :
                        fprintf(stderr,"Unexpected arguement \n");
                        show_help(argv[0]);
                        return 1;
		}
	}

	/* calculate optic vectors & sanity checks */
	B0 = sub(P1, P0);
	B1 = sub(P2, P0);
	B2 = cross(B0,B1);
	if (magnitude(B2) == 0.0) {  //points are colinear
		fprintf(stderr,"Error points are colinear!\n");
		return 1;
	}
	if (dot(B0,B1) != 0.0) { //points are not perpendicular
		fprintf(stderr,"Error points are not perpendicular!\n");
		return 1;
	}

	/* check profile file given*/
	if (mirrorFileName == NULL) {
		fprintf(stderr, "No optical CVS file has been inputed!\n");
		show_help(argv[0]);
		return 1;
	}

	/* check the rows and columns provided */
	if (rows < 1) {
		fprintf(stderr, "Number of rows needs be at least 1!\n");
		show_help(argv[0]);
		return 1;
	}
	if (columns < 1) {
		fprintf(stderr, "Number of columns needs be at least 1!\n");
		show_help(argv[0]);
		return 1;
	}

	/* Load Mirror */
	mZ   = calloc(rows, sizeof(double*));
	for (g = 0; g < rows; g++) {
		mZ[g] = (double *) calloc(columns, sizeof(double));
	} 

	g = read_mirror(mZ, mirrorFileName, rows, columns);
	if (g != 0) return 1; //Error reading profile
	mirrorMax = max_mirror(mZ, rows, columns);
	mirrorMin = min_mirror(mZ, rows, columns);

	/* allocate memory & compute values for lin interp*/
	b00 = calloc(rows, sizeof(double*));
	b01 = calloc(rows, sizeof(double*));
	b10 = calloc(rows, sizeof(double*));
	b11 = calloc(rows, sizeof(double*));
	for (g = 0; g < rows; g++) {
		b00[g] = (double *) calloc(columns, sizeof(double));
		b01[g] = (double *) calloc(columns, sizeof(double));
		b10[g] = (double *) calloc(columns, sizeof(double));
		b11[g] = (double *) calloc(columns, sizeof(double));
	} 
	calc_b00(b00, mZ, rows, columns);
	calc_b01(b01, mZ, rows, columns);
	calc_b10(b10, mZ, rows, columns);
	calc_b11(b11, mZ, rows, columns);

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

	/* Required math */
	B0Unit = unit(B0); //make unit vectors
	B1Unit = unit(B1);
	B2Unit = unit(B2);
	P3 = add(P0,mult(B2Unit,mirrorMax)); //Vector to define "top" of Rectangular cuboid
	P4 = add(P0,mult(B2Unit,mirrorMin)); //Vector to define "bottom" of Rectangular cuboid
	B2 = mult(B2Unit, mirrorMax-mirrorMin); //make lenght of B2 equal to size of Rectangular cuboid
	/* loop over all optical Rays */
	while (!feof(in)) {
		double prop, j, k, theta, VdotN;
		double d0, d1, d2, d3, d4, d5; //distances to walls of Rectangular cuboid plane
		double lmin = 0.0;
		double lmax = 0.0;
		Ray myRay;
		Vector RayMB; //Ray in mirror basis

		/* Read in Ray */
		g = read_ray(in, &myRay);
		if (g != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}

		/* check if 0 intensity skip */
		if (useZero == 0 && myRay.i == 0.0) {
			write_ray(out, myRay);
			continue;
		}

		/* Approximate the optic as a Rectangular cuboid & find the entrence and exit faces*/
		myRay.v = unit(myRay.v); //check that direction is a unit vector
		//North plane
		if (dot(myRay.v, B0) != 0.0) {
			double d, ymin, ymax, yvalue, xmin, xmax, xvalue;
			Vector Q, Vec0, Vec1, Vec2;
			Vec0 = add(P1,P4);
			Vec1 = add(P1,P3);
			Vec2 = add(Vec0,P2);
			d = dot(sub(Vec0, myRay.pos), B0) / dot(myRay.v, B0);
			Q = add(myRay.pos, mult(myRay.v, d)); //location in north plane
			xmin   = dot(B2, Vec0);
			xmax   = dot(B2, Vec1);
			xvalue = dot(B2, Q);
			ymin   = dot(B1, Vec0);
			ymax   = dot(B1, Vec2);
			yvalue = dot(B1, Q);
			if ((xvalue >= xmin) && (xvalue <= xmax) && (yvalue >= ymin) && (yvalue <= ymax)) {
				// point in plane in on Rectangular cuboid face
				if (lmin == lmax) {
					lmin = d;
				}
				else {
					lmax = d;
				}
			}
		}
		//South plane
		if (dot(myRay.v, B0) != 0.0) {
			double d, ymin, ymax, yvalue, xmin, xmax, xvalue;
			Vector Q, Vec0, Vec1, Vec2;
			Vec0 = add(P0,P4);
			Vec1 = add(P0,P3);
			Vec2 = add(Vec0,P2);
			d = dot(sub(Vec0, myRay.pos), B0) / dot(myRay.v, B0);
			Q = add(myRay.pos, mult(myRay.v, d)); //location in north plane
			xmin   = dot(B2, Vec0);
			xmax   = dot(B2, Vec1);
			xvalue = dot(B2, Q);
			ymin   = dot(B1, Vec0);
			ymax   = dot(B1, Vec2);
			yvalue = dot(B1, Q);
			if ((xvalue >= xmin) && (xvalue <= xmax) && (yvalue >= ymin) && (yvalue <= ymax)) {
				// point in plane in on Rectangular cuboid face
				if (lmin == lmax) {
					lmin = d;
				}
				else {
					lmax = d;
				}
			}
		}
		//East plane
		if (dot(myRay.v, B1) != 0.0) {
			double d, ymin, ymax, yvalue, xmin, xmax, xvalue;
			Vector Q, Vec0, Vec1, Vec2;
			Vec0 = add(P0,P4);
			Vec1 = add(P1,P4);
			Vec2 = add(P0,P3);
			d = dot(sub(Vec0, myRay.pos), B1) / dot(myRay.v, B1);
			Q = add(myRay.pos, mult(myRay.v, d)); //location in north plane
			xmin   = dot(B0, Vec0);
			xmax   = dot(B0, Vec1);
			xvalue = dot(B0, Q);
			ymin   = dot(B2, Vec0);
			ymax   = dot(B2, Vec2);
			yvalue = dot(B2, Q);
			if ((xvalue >= xmin) && (xvalue <= xmax) && (yvalue >= ymin) && (yvalue <= ymax)) {
				// point in plane in on Rectangular cuboid face
				if (lmin == lmax) {
					lmin = d;
				}
				else {
					lmax = d;
				}
			}
		}
		//West plane
		if (dot(myRay.v, B1) != 0.0) {
			double d, ymin, ymax, yvalue, xmin, xmax, xvalue;
			Vector Q, Vec0, Vec1, Vec2;
			Vec0 = add(add(P0,P4),P2);
			Vec1 = add(add(P1,P4),P2);
			Vec2 = add(add(P0,P3),P2);
			d = dot(sub(Vec0, myRay.pos), B1) / dot(myRay.v, B1);
			Q = add(myRay.pos, mult(myRay.v, d)); //location in north plane
			xmin   = dot(B1, Vec0);
			xmax   = dot(B1, Vec1);
			xvalue = dot(B1, Q);
			ymin   = dot(B2, Vec0);
			ymax   = dot(B2, Vec2);
			yvalue = dot(B2, Q);
			if ((xvalue >= xmin) && (xvalue <= xmax) && (yvalue >= ymin) && (yvalue <= ymax)) {
				// point in plane in on Rectangular cuboid face
				if (lmin == lmax) {
					lmin = d;
				}
				else {
					lmax = d;
				}
			}
		}
		//Top plane
		if (dot(myRay.v, B2) != 0.0) {
			double d, ymin, ymax, yvalue, xmin, xmax, xvalue;
			Vector Q, Vec0, Vec1, Vec2;
			Vec0 = add(P0,P3);
			Vec1 = add(Vec0,P2);
			Vec2 = add(Vec0,P1);
			d = dot(sub(Vec0, myRay.pos), B1) / dot(myRay.v, B1);
			Q = add(myRay.pos, mult(myRay.v, d)); //location in north plane
			xmin   = dot(B1, Vec0);
			xmax   = dot(B1, Vec1);
			xvalue = dot(B1, Q);
			ymin   = dot(B0, Vec0);
			ymax   = dot(B0, Vec2);
			yvalue = dot(B0, Q);
			if ((xvalue >= xmin) && (xvalue <= xmax) && (yvalue >= ymin) && (yvalue <= ymax)) {
				// point in plane in on Rectangular cuboid face
				if (lmin == lmax) {
					lmin = d;
				}
				else {
					lmax = d;
				}
			}
		}
		//Bottom plane
		if (dot(myRay.v, B2) != 0.0) {
			double d, ymin, ymax, yvalue, xmin, xmax, xvalue;
			Vector Q, Vec0, Vec1, Vec2;
			Vec0 = add(P0,P4);
			Vec1 = add(Vec0,P2);
			Vec2 = add(Vec0,P4);
			d = dot(sub(Vec0, myRay.pos), B2) / dot(myRay.v, B2);
			Q = add(myRay.pos, mult(myRay.v, d)); //location in north plane
			xmin   = dot(B1, Vec0);
			xmax   = dot(B1, Vec1);
			xvalue = dot(B1, Q);
			ymin   = dot(B0, Vec0);
			ymax   = dot(B0, Vec2);
			yvalue = dot(B0, Q);
			if ((xvalue >= xmin) && (xvalue <= xmax) && (yvalue >= ymin) && (yvalue <= ymax)) {
				// point in plane in on Rectangular cuboid face
				if (lmin == lmax) {
					lmin = d;
				}
				else {
					lmax = d;
				}
			}
		}

		RayMB.x = dot(myRay.v,B0);
		RayMB.y = dot(myRay.v,B1);
		RayMB.z = dot(myRay.v,B2);

		//propagate to surface
//		prop = find_surface(myRay, a, b, c, d, N);

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
//		Normal.x = function_eval(dxa, dxb, c, d, N, myRay.pos);
//		Normal.y = function_eval(dya, b, dyc, d, N, myRay.pos);
//		Normal.z = function_eval(dza, b, c, dzd, N, myRay.pos);
//		Normal = unit(Normal);
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
	return 0;
}

