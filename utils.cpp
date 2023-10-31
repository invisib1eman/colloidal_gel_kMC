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
    double alpha=acos(2.0*b-1.0)*step;
    double beta=2.0*M_PI*(c-0.5)*step;
    quarternion rotate=angle_to_quarternion(theta,alpha,beta);
    quarternion neworientation=quartermulti(rotate,old);
    neworientation.normalize();
    return neworientation;
}
