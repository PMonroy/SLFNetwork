#include "vectorXYZ.h"

vectorXYZ::vectorXYZ(void) {  // Void Constructor
	x=0.0;
	y=0.0;
	z=0.0;
}
vectorXYZ::vectorXYZ(double xi,double yi,double zi) {  // Constructor
	x=xi;
	y=yi;
	z=zi;
}
vectorXYZ::vectorXYZ(const vectorXYZ &a){   // Copy Vector Constructor
	x=a.x;
	y=a.y;
	z=a.z;
}

void vectorXYZ::operator=(vectorXYZ a){  // Copy Vector
	x=a.x;
	y=a.y;
	z=a.z;
}
vectorXYZ &vectorXYZ::operator+=(vectorXYZ a) 
{
    x += a.x;
    y += a.y;
    z += a.z;
    return *this;
}
vectorXYZ &vectorXYZ::operator-=(vectorXYZ a) 
{
    x -= a.x;
    y -= a.y;
    z -= a.z;
    return *this;
}
vectorXYZ &vectorXYZ::operator*=(double scalar) 
{
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}
vectorXYZ &vectorXYZ::operator*=(vectorXYZ a) 
{
    x *= a.x;
    y *= a.y;
    z *= a.z;
    return *this;
}
vectorXYZ &vectorXYZ::operator/=(vectorXYZ a) 
{
    x /= a.x;
    y /= a.y;
    z /= a.z;
    return *this;
}


vectorXYZ  vectorXYZ::operator+(vectorXYZ a){  // Addition
  return vectorXYZ(*this)+=a;
}
vectorXYZ  vectorXYZ::operator-(vectorXYZ a){  // Subtraction
  return vectorXYZ(*this)-=a;
}
vectorXYZ  vectorXYZ::operator*(vectorXYZ a){  // Product
  return vectorXYZ(*this)*=a;
}
vectorXYZ  vectorXYZ::operator/(vectorXYZ a){  // Division
  return vectorXYZ(*this)/=a;
}


vectorXYZ::~vectorXYZ(){}    // Destructor


vectorXYZ operator*( double scalar,vectorXYZ a) {
  return vectorXYZ(a) *= scalar;
} 
vectorXYZ operator*(vectorXYZ a, double scalar) {
  return vectorXYZ(a) *= scalar;
} 

ostream &operator<<(ostream &out, vectorXYZ a)     //output
{
        out<<a.x<<" "<<a.y<<" "<<a.z;
        return out;
}

bool operator==(vectorXYZ a,vectorXYZ b)
{
  return (a.x==b.x && a.y==b.y && a.z==b.z);
}
