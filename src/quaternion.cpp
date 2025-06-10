#include "quaternion.h"
quaternion angle_to_quaternion(double theta,double alpha,double beta)
{
  double w=cos(theta*0.5);
  double x=sin(theta*0.5)*sin(alpha)*cos(beta);
  double y=sin(theta*0.5)*sin(alpha)*sin(beta);
  double z=sin(theta*0.5)*cos(alpha);
  return quaternion(w,x,y,z);
}
quaternion quartermulti(quaternion a,quaternion b)
{
  double w_1=a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z;
  double x_1=a.w*b.x+a.x*b.w+(a.y*b.z-a.z*b.y);
  double y_1=a.w*b.y+a.y*b.w+(a.z*b.x-a.x*b.z);
  double z_1=a.w*b.z+a.z*b.w+(a.x*b.y-a.y*b.x);
  return quaternion(w_1,x_1,y_1,z_1);

}
quaternion quartercc(quaternion a)
{
  return quaternion(a.w,-a.x,-a.y,-a.z);
}
quaternion vector_to_quaternion(XYZ a)
{
  return quaternion(0,a.x,a.y,a.z);
}
XYZ quarterrotation(XYZ old,quaternion q)
{
  quaternion quarterv=vector_to_quaternion(old);
  quaternion newquarterv=quartermulti(quartermulti(q,quarterv),quartercc(q));
  return XYZ(newquarterv.x,newquarterv.y,newquarterv.z);
}
quaternion RandomRotate(quaternion old, double step,double a,double b)
{
    
    double theta=step;
    double alpha=acos(2.0*a-1.0);
    double beta=2.0*M_PI*(b-0.5);
    quaternion rotate=angle_to_quaternion(theta,alpha,beta);
    quaternion neworientation=quartermulti(rotate,old);
    neworientation.normalize();
    return neworientation;
}
quaternion RandomRotatestep(double step,double a,double b)
{
    double theta=step;
    double alpha=acos(2.0*a-1.0);
    double beta=2.0*M_PI*(b-0.5);
    quaternion rotate=angle_to_quaternion(theta,alpha,beta);
    return rotate;
}