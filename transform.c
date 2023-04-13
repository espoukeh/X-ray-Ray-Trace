/*
* Transform
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include "ray.h"

static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Transforms the corrdinate system of the incident rays via rotation and\n"
"translation. The rotation is specified by two of the three new axis vectors.\n"
"Note: it first rotates the corrdiate system and then translates the origin.\n"
"If you want to shift then rotate you need to run the program twice.\n"
"Note: all values are metric.\n"
"\n"
"  -h, --help             Display this help message.\n"
"  -i, --input=<file>     Input filename. Default: stdin.\n"
"  -o, --output=<file>    Output filename. Default: stdout.\n"
"      --x='[x,y,z]'      New x vector direction.\n"
"      --y='[x,y,z]'      New y vector direction.\n"
"      --z='[x,y,z]'      New z vector direction.\n"
"                             Note: Only two of x,y,z are required for\n"
"                             a right handed coordiate system.\n"
"      --0='[x,y,z]'      New origin location (default = [0,0,0]).\n" 
"\n");
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
	int c;
	// logic sanity checks
	int x_flag = 0;
	int y_flag = 0;
	int z_flag = 0;
	int o_flag = 0;
	int rotate_flag = 0;
	// input vectors
	Vector x = make_vector(1.0, 0.0, 0.0);
	Vector y = make_vector(0.0, 1.0, 0.0);
	Vector z = make_vector(0.0, 0.0, 1.0);
	Vector o = make_vector(0.0, 0.0, 0.0);
	// math stuff (transform matrix)
	Vector V0, V1, V2; //Transformation matrix as 3 vectors

	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"x",                 1, NULL,               'x'},
		{"y",                 1, NULL,               'y'},
		{"z",                 1, NULL,               'z'},
		{"0",                 1, NULL,               '0'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:x:y:z:0:",
                longopts, NULL)) != -1)
        {

                switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'i' :
			inFileName = strdup(optarg);
			break;

			case 'o' :
			outFileName = strdup(optarg);
			break;

			case 'x' :
			sscanf(optarg,"[%le,%le,%le]", &x.x, &x.y, &x.z);
			x_flag = 1;
			break;

			case 'y' :
			sscanf(optarg,"[%le,%le,%le]", &y.x, &y.y, &y.z);
			y_flag = 1;			
			break;

			case 'z' :
			sscanf(optarg,"[%le,%le,%le]", &z.x, &z.y, &z.z);
			z_flag = 1;
			break;

			case '0' :
			sscanf(optarg,"[%le,%le,%le]", &o.x, &o.y, &o.z);
			o_flag = 1;
			break;

			default :
                        fprintf(stderr,"Unexpected arguement \n");
                        show_help(argv[0]);
                        return 1;

		}
	}
	
	//calculate cross product sanity checks
	if (magnitude(x) == 0.0) {
		fprintf(stderr,"x vector is 0 in magnitude!\n");
	}
	if (magnitude(y) == 0.0) {
		fprintf(stderr,"y vector is 0 in magnitude!\n");
	}
	if (magnitude(z) == 0.0) {
		fprintf(stderr,"z vector is 0 in magnitude!\n");
	}
	if (x_flag + y_flag + z_flag < 0) {
		fprintf(stderr,"unknown error! Blame the programmer!\n");
		return 1;
	}
	if (x_flag + y_flag + z_flag > 3) {
		fprintf(stderr,"unknown error! Blame the programmer!\n");
		return 1;
	}	
	if (x_flag + y_flag + z_flag == 1) {
		fprintf(stderr,"Only 1 roation vector is defined! Please add another direction\n");
		return 1;
	}
	if (x_flag + y_flag + z_flag == 2) {
		if (x_flag == 1) {
			if (y_flag == 1) {
				if (dot(x,y) != 0.0) {
					fprintf(stderr,"Error: x & y basis vector not perpendicular.\n");
					return 1;
				}
				z = cross(x,y);
			} else {
				if (dot(x,z) != 0.0) {
					fprintf(stderr,"Error: x & z basis vector not perpendicular.\n");
					return 1;
				}
				y = cross(z,x);
			}
		} else {
			if (dot(y,z) != 0.0) {
				fprintf(stderr,"Error: y & z basis vector not perpendicular.\n");
				return 1;
			}
			x = cross(y,z);
		}
		rotate_flag = 1;
	}
	if (x_flag + y_flag + z_flag == 3) {
		Vector v;		
		if (dot(x,y) != 0.0) {
			fprintf(stderr,"Error: x & y basis vector not perpendicular.\n");
			return 1;
		}
		if (dot(x,z) != 0.0) {
			fprintf(stderr,"Error: x & z basis vector not perpendicular.\n");
			return 1;
		}
		if (dot(y,z) != 0.0) {
			fprintf(stderr,"Error: y & z basis vector not perpendicular.\n");
			return 1;
		}
		v = cross(x,y);
		if ((v.x == z.x) && (v.y == z.y) && (v.z == z.z)) {
			rotate_flag = 1;
		}
		else {
			fprintf(stderr,"Error: you are unsing a left handed corrdinate system.\n");
			return 1;
		}
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

	/* calculate transformation matrix */
	if (rotate_flag == 1) {

		V0 = unit(x);
		V1 = unit(y);
		V2 = unit(z);

		/* For a given vecotr v the following is true:
		x' = dot(V0, v);
		y' = dot(V1, v);
		z' = dot(V2, v);
		*/

//		V0.x = (y.z*z.y - y.y*z.z)/(x.z*y.y*z.x - x.y*y.z*z.x - x.z*y.x*z.y + x.x*y.z*z.y + x.y*y.x*z.z - x.x*y.y*z.z);
//		V0.y = (y.x*z.z - y.z*z.x)/(x.z*y.y*z.x - x.y*y.z*z.x - x.z*y.x*z.y + x.x*y.z*z.y + x.y*y.x*z.z - x.x*y.y*z.z);
//	 	V0.z = (y.y*z.x - y.x*z.y)/(x.z*y.y*z.x - x.y*y.z*z.x - x.z*y.x*z.y + x.x*y.z*z.y + x.y*y.x*z.z - x.x*y.y*z.z);
//		V0 = unit(V0);
	}

	/* loop over vectors */
	while (!feof(in)) {				
		Ray myRay;
		c = read_ray(in, &myRay);
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}
		
		/* check if need to rotate */
		if (rotate_flag == 1) {
			Vector new_pos, new_v;
			new_pos.x = dot(V0, myRay.pos);
			new_pos.y = dot(V1, myRay.pos);
			new_pos.z = dot(V2, myRay.pos);
			new_v.x = dot(V0, myRay.v);
			new_v.y = dot(V1, myRay.v);
			new_v.z = dot(V2, myRay.v);
			myRay.pos = new_pos;
			myRay.v = new_v;
		}

		/* check if need to translate */
		if (o_flag == 1) {
			myRay.pos = sub(myRay.pos, o);
		}

		//print ray
		write_ray(out, myRay);
	}

	fclose(in);
	fclose(out);
	return 0;
}
