/*
* grep_rays
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


enum {
	EQ,
	NE,
	LE,
	LT,
	GT,
	GE
};

static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Filters ray data to satisfy specific criteria in one of the data fields.\n"
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
"  -l, --lt=<num>          Seletion criteria is values less than <num>.\n"
"  -a, --le=<num>          Seletion criteria is values less than or equal to <num>.\n"
"  -g, --gt=<num>          Seletion criteria is values greater than <num>.\n"
"  -b, --ge=<num>          Seletion criteria is values greater than or equal to <num>.\n"
"  -e, --eq=<num>          Seletion criteria is values equal to <num>.\n"
"  -n, --ne=<num>          Seletion criteria is values not equal to <num>.\n"
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
	int criteria;
	int sanityCheck = 0;
	unsigned long  criteriaUL;
	double criteriaDouble;
	
	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"field",             1, NULL,               'f'},
		{"lt",                1, NULL,               'l'},
		{"le",                1, NULL,               'a'},
		{"gt",                1, NULL,               'g'},
		{"ge",                1, NULL,               'b'},
		{"eq",                1, NULL,               'e'},
		{"ne",                1, NULL,               'n'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:a:b:e:f:g:l:n:",
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

			case 'l' :
			valueString = strdup(optarg);
			criteria = LT;
			sanityCheck++;
			break;			

			case 'a' :
			valueString = strdup(optarg);
			criteria = LE;
			sanityCheck++;
			break;

			case 'g' :
			valueString = strdup(optarg);
			criteria = GT;
			sanityCheck++;
			break;

			case 'b' :
			valueString = strdup(optarg);
			criteria = GE;
			sanityCheck++;
			break;

			case 'e' :
			valueString = strdup(optarg);
			criteria = EQ;
			sanityCheck++;
			break;

			case 'n' :
			valueString = strdup(optarg);
			criteria = NE;
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
		fprintf(stderr,"Error only one selection criteria can be chosen\n");
		return 1;
	}
	if (field == -1) {
		fprintf(stderr,"Error a field has not been selected\n");
		return 1;
	}

	/* convert value to apprpriate form*/
	if (field == RN) {
		criteriaUL = atoi(valueString);
		
	} else {
		criteriaDouble = atof(valueString);
	}

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
		double valueDouble;
		unsigned long valueUL;
		c = read_ray(in, &myRay);
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}

		switch (field) {
			case RN:
			valueUL = myRay.ray;
			break;

			case PX:
			valueDouble = myRay.pos.x;
			break;
			
			case PY:
			valueDouble = myRay.pos.y;
			break;	

			case PZ:
			valueDouble = myRay.pos.z;
			break;	

			case DX:
			valueDouble = myRay.v.x;
			break;
			
			case DY:
			valueDouble = myRay.v.y;
			break;	

			case DZ:
			valueDouble = myRay.v.z;
			break;	

			case W:
			valueDouble = myRay.w;
			break;
			
			case PL:
			valueDouble = myRay.l;
			break;	

			case I:
			valueDouble = myRay.i;
			break;	

			case PO:
			valueDouble = myRay.p;
			break;
		}

 		/* not positive about casting unsigned long to double... 
		so just being safe (or stupid) */
		if (field == RN) { 
			switch (criteria) {
				case EQ:
				if (valueUL == criteriaUL) {
					write_ray(out, myRay);
				}
				break;
				
				case NE:
				if (valueUL != criteriaUL) {
					write_ray(out, myRay);
				}
				break;

				case LT:
				if (valueUL < criteriaUL) {
					write_ray(out, myRay);
				}
				break;
				
				case LE:
				if (valueUL <= criteriaUL) {
					write_ray(out, myRay);
				}
				break;
				
				case GT:
				if (valueUL > criteriaUL) {
					write_ray(out, myRay);
				}
				break;
				
				case GE:
				if (valueUL >= criteriaUL) {
					write_ray(out, myRay);
				}
				break;
			}
		}
		else {
			switch (criteria) {
				case EQ:
				if (valueDouble == criteriaDouble) {
					write_ray(out, myRay);
				}
				break;
				
				case NE:
				if (valueDouble != criteriaDouble) {
					write_ray(out, myRay);
				}
				break;

				case LT:
				if (valueDouble < criteriaDouble) {
					write_ray(out, myRay);
				}
				break;
				
				case LE:
				if (valueDouble <= criteriaDouble) {
					write_ray(out, myRay);
				}
				break;
				
				case GT:
				if (valueDouble > criteriaDouble) {
					write_ray(out, myRay);
				}
				break;
				
				case GE:
				if (valueDouble >= criteriaDouble) {
					write_ray(out, myRay);
				}
				break;
			}
		}	
	}
	fclose(in);
	fclose(out);
	return 0;
}
