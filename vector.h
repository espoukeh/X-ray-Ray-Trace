/*
* Math needed for simple 3D vector analysis
*/

typedef struct {
	double x;
	double y;
	double z;
} Vector; 

Vector make_vector(double x, double y, double z);
Vector add(Vector u, Vector v);
Vector mult(Vector v, double s);
Vector sub(Vector u, Vector v);
Vector cross(Vector u, Vector v);
Vector unit(Vector v);
double dot(Vector u, Vector v);
double magnitude(Vector v);

