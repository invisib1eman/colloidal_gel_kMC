//NANOROD: utils.cpp Utilities Functions (Revision Date: Oct 27, 2023)
#include "utils.h"
//do not assume where the particle is (particle may have crossed the box twice or more)
XYZ image(XYZ p, double L)
{
  XYZ xyz;
  xyz.x=myfmod(p.x+0.5*L, L)-0.5*L;
  xyz.y=myfmod(p.y+0.5*L, L)-0.5*L;
  xyz.z=myfmod(p.z+0.5*L, L)-0.5*L;
  return xyz;
}
double myfmod(double x, double y)
{
  double temp=fmod(x,y);
  if(temp>=0.0)
    return temp;
  else
    return temp+y;
}
//Minimum image distance squared: pass original coordinates
double min_d2(XYZ a, XYZ b, double L)
{
    a=image(a,L);
    b=image(b,L);
    XYZ d=a-b;
    d.my_abs();
    if(d.x>=0.5*L)
      d.x=L-d.x;
    
    if(d.y>=0.5*L)
      d.y=L-d.y;
    
    if(d.z>=0.5*L)
      d.z=L-d.z;
    
    return d.norm2();
}
XYZ RandomTranslate(XYZ old, double step,double u,double v)
{
    
    double theta=2.0*M_PI*u;
    double phi=acos(2.0*v-1.0);
    return XYZ(old.x+step*cos(theta)*sin(phi),old.y+step*sin(theta)*sin(phi),old.z+step*cos(phi));
}


quarternion RandomRotate(quarternion old, double step,double a,double b,double c)
{
    
    double theta=2.0*M_PI*a*step;
    double alpha=acos(2.0*b-1.0);
    double beta=2.0*M_PI*(c-0.5);
    quarternion rotate=angle_to_quarternion(theta,alpha,beta);
    quarternion neworientation=quartermulti(rotate,old);
    neworientation.normalize();
    return neworientation;
}
int getNum(vector<int>& v)
{
 
    // Size of the vector
    int n = v.size();
 
    // Generate a random number
    srand(time(NULL));
 
    // Make sure the number is within
    // the index range
    int index = rand() % n;
 
    // Get random number from the vector
    int num = v[index];
 
    // Remove the number from the vector
    swap(v[index], v[n - 1]);
    v.pop_back();
 
    // Return the removed number
    return num;
}
 
// Function to generate n non-repeating random numbers
vector<int> generateRandom(int n)
{
    vector<int> v(n);
 
    // Fill the vector with the values
    // 1, 2, 3, ..., n
    for (int i = 0; i < n; i++)
        v[i] = i + 1;
    vector<int> randomv;
    // While vector has elements
    // get a random number from the vector and print it
    while (v.size()) {
        randomv.push_back(getNum(v));
    }
    return randomv;
}
//inner product of two vectors
double inner_product(XYZ a,XYZ b)
{
  return a.x*b.x+a.y*b.y+a.z*b.z;
}
//cross product of two vectors
XYZ cross_product(XYZ a,XYZ b)
{
  return XYZ(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}
//calculate angle between two vectors
double angle_vectors(XYZ a,XYZ b)
{
  return acos(inner_product(a,b)/(a.norm()*b.norm()));
}
double dihedral_vectors(XYZ a,XYZ b,XYZ c)
{
    return atan2(inner_product(cross_product(a,b),cross_product(b,c))/(cross_product(a,b).norm()*cross_product(b,c).norm()),b.norm()*inner_product(a,cross_product(b,c))/(cross_product(a,b).norm()*cross_product(b,c).norm()));
}