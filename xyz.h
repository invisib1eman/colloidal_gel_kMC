//NANOROD: XYZ Class to store and manipulate coordinates  (Revision Date: Oct 27, 2023)
#ifndef _XYZ_H
#define _XYZ_H
#include "header.h"
class XYZ
{
  public:
      double x,y,z;
      void set(double _x, double _y, double _z){x=_x;        y=_y;   z=_z;}
      XYZ(double _x, double _y, double _z){x=_x;        y=_y;   z=_z;}
      XYZ(){x=0.0;      y=0.0;  	z=0.0;}
      XYZ operator + (XYZ& other){return XYZ(x+other.x,y+other.y,z+other.z);}
      XYZ operator - (XYZ& other){return XYZ(x-other.x,y-other.y,z-other.z);}
      XYZ operator / (double s){return XYZ(x/s,y/s,z/s);}
      double real_d2(XYZ& other){return (x-other.x)*(x-other.x)+(y-other.y)*(y-other.y)+(z-other.z)*(z-other.z);}
      double real_d(XYZ& other){return sqrt((x-other.x)*(x-other.x)+(y-other.y)*(y-other.y)+(z-other.z)*(z-other.z));}
      void my_abs(){x=fabs(x); y=fabs(y); z=fabs(z);}
      double norm2(){ return x*x+y*y+z*z; }

};
#endif
