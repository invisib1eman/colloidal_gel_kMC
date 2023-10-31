//https://graphics.stanford.edu/courses/cs348a-17-winter/Papers/quaternion.pdf
//https://theswissbay.ch/pdf/Gentoomen%20Library/Game%20Development/Programming/Graphics%20Gems%203.pdf
#ifndef _quarternion_H
#define _quarternion_H
#include "header.h"
#include "xyz.h"

class quarternion
{
    public:
        double w,x,y,z;
        void set(double _w, double _x, double _y, double _z){w=_w;        x=_x;   y=_y;z=_z;}
        quarternion(double _w,double _x, double _y, double _z){w=_w;        x=_x;   y=_y;z=_z;}
        quarternion(){w=0;        x=0;   y=0;z=0;}
        
        quarternion operator + (quarternion& other){return quarternion(w+other.w,x+other.x,y+other.y,z+other.z);}
        quarternion operator - (quarternion& other){return quarternion(w-other.w,x-other.x,y-other.y,z-other.z);}
        quarternion operator * (double s){return quarternion(w*s,x*s,y*s,z*s);}
        quarternion operator / (double s){return quarternion(w/s,x/s,y/s,z/s);}
        void normalize()
        {double norm=w*w+x*x+y*y+z*z;
         w=w/norm;
         x=x/norm;
         y=y/norm;
         z=z/norm;
        }
        
         
};
quarternion angle_to_quarternion(double theta,double alpha,double beta);
quarternion quartermulti(quarternion a,quarternion b);
quarternion quartercc(quarternion a);
quarternion vector_to_quarternion(XYZ a);
XYZ quarterrotation(XYZ old,quarternion q);
#endif