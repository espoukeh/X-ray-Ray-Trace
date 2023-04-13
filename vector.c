#include <math.h>
#include "vector.h"

Vector make_vector(double x, double y, double z)
{
	Vector v;
	v.x = x;
	v.y = y;
	v.z = z;
	return v;
}

Vector add(Vector u, Vector v) 
{
	Vector w;
	w.x = u.x + v.x;
	w.y = u.y + v.y;
	w.z = u.z + v.z;
	return w;
}

Vector mult(Vector v, double s)
{
	Vector w;
	w.x = v.x * s;
	w.y = v.y * s;
	w.z = v.z * s;
	return w;
}

//subtraction is syntatic sugar, but might as well include it
Vector sub(Vector u, Vector v) 
{
	Vector w;
	w.x = u.x - v.x;
	w.y = u.y - v.y;
	w.z = u.z - v.z;
	return w;
}

Vector cross(Vector u, Vector v)
{
	Vector w;
	w.x = (u.y * v.z) - (u.z * v.y);
	w.y = (u.z * v.x) - (u.x * v.z);
	w.z = (u.x * v.y) - (u.y * v.x);
	return w;
}

double dot(Vector u, Vector v)
{
	double s;
	s = (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
	return s;
}

double magnitude(Vector v)
{
	double s;
	s = dot(v,v);
	s = sqrt(s);
	return s;
}

Vector unit(Vector v)
{
	return mult(v, 1.0/magnitude(v));
}

