//NANOROD: utils.cpp Utilities Functions (Revision Date: Oct 27, 2023)
#include "utils.h"
//do not assume where the particle is (particle may have crossed the box twice or more)
XYZ image(XYZ p, double BoxLength)//note that every time we calculate vectors or distances we need to image XYZ
{
  XYZ xyz;
  xyz.x=myfmod(p.x+0.5*BoxLength, BoxLength)-0.5*BoxLength;
  xyz.y=myfmod(p.y+0.5*BoxLength, BoxLength)-0.5*BoxLength;
  xyz.z=myfmod(p.z+0.5*BoxLength, BoxLength)-0.5*BoxLength;
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
double min_d2(XYZ a, XYZ b, double BoxLength)
{
    a=image(a,BoxLength);
    b=image(b,BoxLength);
    XYZ d=a-b;
    d.my_abs();
    if(d.x>=0.5*BoxLength)
      d.x=BoxLength-d.x;
    
    if(d.y>=0.5*BoxLength)
      d.y=BoxLength-d.y;
    
    if(d.z>=0.5*BoxLength)
      d.z=BoxLength-d.z;
    
    return d.norm2();
}
XYZ real_vector(XYZ origin,double BoxLength)//after image,get the real vectors
{
  XYZ real=origin;
  if(origin.x>0.5*BoxLength)
  {
    real.x=origin.x-BoxLength;
  }
  if(origin.x<-0.5*BoxLength)
  {
    real.x=origin.x+BoxLength;
  }
  if(origin.y>0.5*BoxLength)
  {
    real.y=origin.y-BoxLength;
  }
  if(origin.y<-0.5*BoxLength)
  {
    real.y=origin.y+BoxLength;
  }
  if(origin.z>0.5*BoxLength)
  {
    real.z=origin.z-BoxLength;
  }
  if(origin.z<-0.5*BoxLength)
  {
    real.z=origin.z+BoxLength;
  }
  return real;
}

XYZ RandomTranslate(XYZ old, double step,double u,double v)
{
    
    double theta=2.0*M_PI*u;
    double phi=acos(2.0*v-1.0);
    return XYZ(old.x+step*cos(theta)*sin(phi),old.y+step*sin(theta)*sin(phi),old.z+step*cos(phi));
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
    return atan2(b.norm()*inner_product(a,cross_product(b,c))/(cross_product(a,b).norm()*cross_product(b,c).norm()),inner_product(cross_product(a,b),cross_product(b,c))/(cross_product(a,b).norm()*cross_product(b,c).norm()));
}
int GridIndex_index(int i,int j,int k,int n)
{
  if(i>n-1)
    i-=n;
  else if(i<0)
    i+=n;
  if(j>n-1)
    j-=n;
  else if(j<0)
    j+=n;
  if(k>n-1)
    k-=n;
  else if(k<0)
    k+=n;
  return n*n*k+n*j+i;
}
int GridIndex_xyz(XYZ& p,int n,double dl,double BoxLength)
{
  int i=int(floor((p.x+0.5*BoxLength)/dl));
  int j=int(floor((p.y+0.5*BoxLength)/dl));
  int k=int(floor((p.z+0.5*BoxLength)/dl));
  return n*n*k+n*j+i;

}
void GridLoc(int& i,int& j,int& k,int n,int index)
{
  i=index%n;
  j=((index-i)/n)%n;
  k=(index-i-j*n)/(n*n);
}
int neighborarm(int n)
{
  return n+1-2*(n%2);
}
// vector<Molecule> generate_unitcell(Molecule M1)
// {

// }
// vector<Molecule> generate_newunit(vector<Molecule> origin_cell,int latticeindex1,int latticeindex2,int latticeindex3)
// {
  
// }
