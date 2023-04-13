/*
* optic_baffle
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
"Creates an optic baffle plane used to define edges of the optic or walls.\n"
"Sets all rays outside boundry to 0 intensity. Baffles are repersented by a\n"
"plane where the normal points towards the propagation side:\n"
"  n=(p1-p0) cross (p2-p0).\n"
"Notes: all values are metric.\n"
"       to swithc acceptable side swap p1 and p2 values\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"      --P0='[x,y,z]'      First point (default = [0,0,0].\n"
"      --P1='[x,y,z]'      Second point (default = [1,0,0].\n"
"      --p2='[x,y,z]'      Third point (default = [0,1,0].\n"
"\n");
}

int main(int argc, char *argv[])
{

	/* initial variable */
	FILE *in = stdin;
	FILE *out = stdout;
	char *inFileName = NULL;
	char *outFileName = NULL;
	Vector P0 = make_vector(0.0, 0.0, 0.0);
	Vector P1 = make_vector(1.0, 0.0, 0.0);
	Vector P2 = make_vector(0.0, 1.0, 0.0);
	Vector N; //normal vector
	double value, len;
	int c;

/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"P0",                1, NULL,               'a'},
		{"P1",                1, NULL,               'b'},
		{"P2",                1, NULL,               'c'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:a:b:c:",
                longopts, NULL)) != -1)
        {

                switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'i' :
			inFileName = strdup(optarg);
			break;

			case 'f' :
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

                        default :
                        fprintf(stderr,"Unexpected arguement \n");
                        show_help(argv[0]);
                        return 1;

		}
	}
	
	//calculate cross product
	N = cross(sub(P1,P0),sub(P2,P0));
	if (magnitude(N) == 0.0) {  //see if points are colinear
		fprintf(stderr,"Error points are colinear!\n");
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

	/* loop over vectors */
	while (!feof(in)) {
		Ray myRay;
		c = read_ray(in, &myRay);
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}
		myRay.v = unit(myRay.v); //check that direction is a unit vector
		if (dot(myRay.v,N) == 0.0) {
			fprintf(stderr,"Warning: Ray %lu propagates parallel to"
			 " plane!\n Skipping ray, but copying it to output.\n", myRay.ray);
		} else {	
			//caluclate distance from the location point in myRay to 
			//plane along propagation direction
			len = dot(sub(P0,myRay.pos),N)/dot(myRay.v,N);
			//calculate which side of plane ray is on
			value = dot(N,mult(myRay.v,len));
			if (value >= 0.0) { //Ray on wrong side of plane
				myRay.i = 0.0; //set intensity to 0
			}
		}	
		//print ray
		write_ray(out, myRay);
	}
	fclose(in);
	fclose(out);
	return 0;
}
	
