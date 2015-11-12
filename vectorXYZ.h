#ifndef VECTORXYZ
#define VECTORXYZ

#include <iostream>

using namespace std;

class  vectorXYZ 
{

 public:
  double x,y,z; // Order list of 3 elements |x|y|z|

  vectorXYZ(void);                          // Zero Vector Constructor
  vectorXYZ(double xi,double yi,double zi); // Constructor
  vectorXYZ(const vectorXYZ &a);            // Copy Vector Constructor


  void operator=(vectorXYZ a);
  vectorXYZ &operator+=(vectorXYZ a); 
  vectorXYZ &operator-=(vectorXYZ a); 
  vectorXYZ &operator*=(vectorXYZ a);
  vectorXYZ &operator*=(double scalar); 
  vectorXYZ &operator/=(vectorXYZ a);

  vectorXYZ operator+( vectorXYZ a);      // Addition
  vectorXYZ operator-( vectorXYZ a);      // Subtraction
  vectorXYZ operator*( vectorXYZ a);      // Product element by element
  vectorXYZ operator/( vectorXYZ a);      // Product element by element

  friend ostream &operator<<(ostream &out, vectorXYZ a);

  ~vectorXYZ();	                            // Destructor
};

vectorXYZ  operator*(double scalar, vectorXYZ a); // Product element by element
vectorXYZ  operator*(vectorXYZ a,double scalar);  // Product element by element

bool operator==(vectorXYZ a,vectorXYZ b);

#endif 

/*
  struct  vectorXYZ {
  double x;
  double y;
  double z;
  };
  
  EXAMPLE OF ALLOCATION
  
  int i;
  p1 = new Vector3D * [10];
  for(i=0; i< 10; i++)
  p1[i] = new Vector3D [5];
  
*/

