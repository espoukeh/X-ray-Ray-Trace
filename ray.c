#include <stdio.h>
#include "ray.h"

Ray read_ray_string(char *line)
{
	Ray r;
	double x,y,z,vx,vy,vz,w,l,i,p;
	long unsigned int ray;
	sscanf(line, "%lu, %le, %le, %le, %le, %le, %le, %le, %le, %le, %le \n", 
	  &ray, &x, &y, &z, &vx, &vy, &vz, &w, &l, &i, &p);
	r.pos = make_vector( x, y, z );
	r.v = make_vector( vx, vy, vz );
	r.w = w;
	r.i = i;
	r.l = l;
	r.p = p;
	r.ray = ray;
	return r;
}

void write_ray_string(char *line, Ray r, int p )
{
	char format[RAY_STRING_LENGTH];	
	sprintf(format,"%%lu, %%1.%ile, %%1.%ile, %%1.%ile, %%1.%ile, %%1.%ile, "
	   "%%1.%ile, %%1.%ile, %%1.%ile, %%1.%ile, %%1.%ile \n", p, p, p, p, p, 
	   p, p, p, p, p);
	sprintf(line, format, r.ray, r.pos.x, r.pos.y, r.pos.z, r.v.x, r.v.y,
	   r.v.z, r.w, r.l, r.i, r.p);
}

void write_header(char *line)
{
	sprintf(line,"ray #, position x, position y, position z, direction x, "
	   "direction y, direction z, wavelength, path lenght, intensity, "
	   "polarization\n");
}

void write_test(char *line, Ray r)
{
	sprintf(line,"%lu, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
	   r.ray, r.pos.x, r.pos.y, r.pos.z, r.v.x, r.v.y, r.v.z, r.w, r.l, r.i, r.p);
}

int write_ray(FILE *fh, Ray r)
{
	double data[] = {r.pos.x, r.pos.y, r.pos.z, r.v.x, r.v.y, r.v.z, r.w,
	                 r.l, r.i, r.p};
	if( fwrite( &(r.ray),   sizeof(unsigned long),  1, fh) !=  1) return 1;
	if( fwrite( data,       sizeof(double),        10, fh) != 10) return 1;
	return 0;
}

int read_ray(FILE *fh, Ray *r)
{
	double data[10];
	unsigned long i;
	if( fread( &i, sizeof(unsigned long),  1, fh) !=  1) return 1;
	if( fread( data,       sizeof(double),        10, fh) != 10) return 1;
	r->ray = i;
	r->pos.x = data[0];
	r->pos.y = data[1];
	r->pos.z = data[2];
	r->v.x   = data[3];
	r->v.y   = data[4];
	r->v.z   = data[5];
	r->w     = data[6];
	r->l     = data[7];
	r->i     = data[8];
	r->p     = data[9];
	return 0;
}
