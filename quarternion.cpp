#include "quarternion.h"
quarternion angle_to_quarternion(double theta,double alpha,double beta)
{
  double w=cos(theta*0.5);
  double x=sin(theta*0.5)*sin(alpha)*cos(beta);
  double y=sin(theta*0.5)*sin(alpha)*sin(beta);
  double z=sin(theta*0.5)*cos(alpha);
  return quarternion(w,x,y,z);
}
quarternion quartermulti(quarternion a,quarternion b)
{
  double w_1=a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z;
  double x_1=a.w*b.x+a.x*b.w+(a.y*b.z-a.z*b.y);
  double y_1=a.w*b.y+a.y*b.w+(a.z*b.x-a.x*b.z);
  double z_1=a.w*b.z+a.z*b.w+(a.x*b.y-a.y*b.x);
  return quarternion(w_1,x_1,y_1,z_1);

}
quarternion quartercc(quarternion a)
{
  return quarternion(a.w,-a.x,-a.y,-a.z);
}
quarternion vector_to_quarternion(XYZ a)
{
  return quarternion(0,a.x,a.y,a.z);
}
XYZ quarterrotation(XYZ old,quarternion q)
{
  quarternion quarterv=vector_to_quarternion(old);
  quarternion newquarterv=quartermulti(quartermulti(q,quarterv),quartercc(q));
  return XYZ(newquarterv.x,newquarterv.y,newquarterv.z);
}

