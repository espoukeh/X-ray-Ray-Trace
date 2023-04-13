/*
* Needed to define a ray interface
*/

#include "vector.h"

//define the max lenght of a string required to hold a Ray
#define RAY_STRING_LENGTH  1024 

//define what goes in an array
typedef struct {
	Vector pos; //position
	Vector v; //direction
	double w; //wavelength
	double l; //pathLength
	double i; //Intensity
	double p; //Polarization
	unsigned long ray; //ray number
} Ray; 

// methods to read and write Rays to make other programs easier
Ray read_ray_string(char *line);
void write_ray_string(char *line, Ray r, int p);
void write_test(char *line, Ray r);
void write_header(char *line);
int write_ray(FILE *fh, Ray r);
int read_ray(FILE *fh, Ray * r);
