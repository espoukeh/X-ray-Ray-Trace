/*
* set_rays
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include "ray.h"

#define RN  1
#define PX  2
#define PY  3
#define PZ  4
#define DX  5
#define DY  6
#define DZ  7
#define W   8
#define PL  9
#define I  10
#define PO 11


static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"sets the value of a data field for all rays to a specified value.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"  -f, --field=<num>       Field identiy number to parce data on:\n"
"                              1 =  Ray number\n"
"                              2 =  Position  x\n"
"                              3 =  Position  y\n"
"                              4 =  Position  z\n"
"                              5 =  Direction x\n"
"                              6 =  Direction y\n"
"                              7 =  Direction z\n"
"                              8 =  Wavelength\n"
"                              9 =  Path length\n"
"                             10 =  Intensity\n"
"                             11 =  Polarization\n"
"  -v, --value=<num>       Value of selected field equails <num>.\n"
"\n");
}

int main(int argc, char *argv[])
{
	/* initial variable */
	FILE *in = stdin;
	FILE *out = stdout;
	char *inFileName = NULL;
	char *outFileName = NULL;
	char *valueString = NULL;
	int c;
	int field = -1;
	int sanityCheck = 0;
	unsigned long  newUL;
	double newDouble;
	
	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"field",             1, NULL,               'f'},
		{"value",             1, NULL,               'v'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:f:v:",
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

			case 'f' :
			field = atoi(optarg);
			if (field < 1 || field >11) {
				fprintf(stderr,"Field value out of range\n");		
				return 1;
			}
			break;

			case 'v' :
			valueString = strdup(optarg);
			sanityCheck++;
			break;			

                        default :
                        fprintf(stderr,"Unexpected arguement \n");
                        show_help(argv[0]);
                        return 1;
		}
	}

	/* sanity check */
	if (sanityCheck != 1) {  
		fprintf(stderr,"Error a value has not been specified\n");
		return 1;
	}
	if (field == -1) {
		fprintf(stderr,"Error a field has not been selected\n");
		return 1;
	}

	/* convert value to apprpriate form*/
	if (field == RN) {
		newUL = atoi(valueString);
		
	} else {
		newDouble = atof(valueString);
	}
        free(valueString);

	/* open file(s) if necessary */
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

	while (!feof(in)) {
		Ray myRay;
		c = read_ray(in, &myRay);
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}

		switch (field) {
			case RN:
			myRay.ray = newUL;
			break;

			case PX:
			myRay.pos.x = newDouble;
			break;
			
			case PY:
			myRay.pos.y = newDouble;
			break;	

			case PZ:
			myRay.pos.z = newDouble;
			break;	

			case DX:
			myRay.v.x = newDouble;
			break;
			
			case DY:
			myRay.v.y = newDouble;
			break;	

			case DZ:
			myRay.v.z = newDouble;
			break;	

			case W:
			myRay.w = newDouble;
			break;
			
			case PL:
			myRay.l = newDouble;
			break;	

			case I:
			myRay.i = newDouble;
			break;	

			case PO:
			myRay.p = newDouble;
			break;
		}

		write_ray(out,myRay);
	}
	fclose(in);
	fclose(out);
	return 0;
}
