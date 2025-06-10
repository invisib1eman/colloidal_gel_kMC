
//https://graphics.stanford.edu/courses/cs348a-17-winter/Papers/quaternion.pdf
//https://theswissbay.ch/pdf/Gentoomen%20Library/Game%20Development/Programming/Graphics%20Gems%203.pdf
#ifndef _quaternion_H
#define _quaternion_H
#include "header.h"
#include "quaternion.h"
#include "xyz.h"

class quaternion
{
    public:
        double w,x,y,z;
        void set(double _w, double _x, double _y, double _z){w=_w;        x=_x;   y=_y;z=_z;}
        quaternion(double _w,double _x, double _y, double _z){w=_w;        x=_x;   y=_y;z=_z;}
        quaternion(){w=0;        x=0;   y=0;z=0;}
        quaternion operator + (quaternion& other){return quaternion(w+other.w,x+other.x,y+other.y,z+other.z);}
        quaternion operator - (quaternion& other){return quaternion(w-other.w,x-other.x,y-other.y,z-other.z);}
        quaternion operator * (double s){return quaternion(w*s,x*s,y*s,z*s);}
        quaternion operator / (double s){return quaternion(w/s,x/s,y/s,z/s);}
        void normalize()
        {double norm=sqrt(w*w+x*x+y*y+z*z);
         w=w/norm;
         x=x/norm;
         y=y/norm;
         z=z/norm;
        }
        
         
};
quaternion angle_to_quaternion(double theta,double alpha,double beta);
quaternion quartermulti(quaternion a,quaternion b);
quaternion quartercc(quaternion a);
quaternion vector_to_quaternion(XYZ a);
XYZ quarterrotation(XYZ old,quaternion q);
quaternion RandomRotate(quaternion old, double step,double a,double b);
quaternion RandomRotatestep(double step,double a,double b);
#endif
