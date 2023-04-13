/*
* Circ_aperture
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
"Creates a circular aperture and propagates all rays to the aperture.\n"
"Sets all rays outside aperture to 0 intensity. Aperture defined by the \n"
"center position C, the radius  R (on the aperture), and the normal\n"
"(perpendicular) vector to the aperture N.\n" 
"Note: all values are metric.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"      --C='[x,y,z]'       Center of aperture (default = [0,0,0]).\n"
"      --R=r               Radius of aperture (default = 1.0).\n"
"      --N='[x,y,z]'       Vector defining normal (default = [0,1,0]).\n"
"      --invert            Inverts the aperture making a beam stop.\n"
"  -z, --zero              Enables propagation of 0 intensity rays\n"
"                            By default not propaged.\n" 
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
	int invert = 0;
	int useZero = 0;
	Vector Center = make_vector(0.0, 0.0, 0.0);
	Vector Normal = make_vector(0.0, 1.0, 0.0);
	double Radius = 1.0;

	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"N",                 1, NULL,               'a'},
		{"R",                 1, NULL,               'b'},
		{"C",                 1, NULL,               'c'},
		{"invert",            0, NULL,               'f'},
		{"zero",              0, NULL,               'z'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:a:b:c:fz",
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

			case 'a' :
			sscanf(optarg,"[%le,%le,%le]", &Normal.x, &Normal.y, &Normal.z);
			break;

			case 'b' :
			Radius = atof(optarg);
			break;

			case 'c' :
			sscanf(optarg,"[%le,%le,%le]", &Center.x, &Center.y, &Center.z);
			break;

			case 'f' :
			invert = 1;
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

	Normal = unit(Normal); // make positive N is a unit normal vector
	
	/* loop over vectors */
	while (!feof(in)) {				
		double distance;
		Ray myRay;
		c = read_ray(in, &myRay);
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}
		myRay.v = unit(myRay.v); //check that direction is a unit vector
		
		/* check if using 0 or skip */
		if (useZero == 0 && myRay.i == 0.0) {
			write_ray(out, myRay);
			continue;
		}

		/* check ray is not parallel to aperture */
		if (dot(myRay.v,Normal) == 0.0) {
			fprintf(stderr,"Warning: Ray %lu propagates parallel to"
			 " aperture!\n Skipping ray, but copying it to output.\n", myRay.ray);
			write_ray(out, myRay);
			continue;
		} 

		//caluclate distance from the location point in myRay to 
		//aperture plane along propagation direction
		distance = dot(sub(Center,myRay.pos),Normal)/dot(myRay.v,Normal);
		if (distance < 0.0) {
			fprintf(stderr,"Warning Ray %lu is furthern then the aperture"
			 " back propagating to aperture\n", myRay.ray);
		}
		myRay.pos = add(myRay.pos,mult(myRay.v,distance));
		myRay.l = myRay.l + distance;

		//calculate where myRay.pos is with respect to aperture center and if inside radius
		if (Radius >= magnitude(sub(myRay.pos,Center))) {
			//in aperture
			if (invert == 1) myRay.i = 0.0; // is a beam block & inside beam block
		} else {
			if (invert == 0) myRay.i = 0.0; // is an aperture & outsid aperture
		}

		//print ray
		write_ray(out, myRay);
	}

	fclose(in);
	fclose(out);
	return 0;
}
