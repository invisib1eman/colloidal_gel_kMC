//NANOROD: hbond Class to manage the bond information  (Revision Date: Oct 27, 2023)
#ifndef _hbond_H
#define _hbond_H
#include "header.h"
class hbond
{
  public:
      int M1,M2,arm1,arm2;
      void set(int _M1, int _M2, int _arm1,int _arm2){M1=_M1;        M2=_M2;   arm1=_arm1;arm2=_arm2;}
      hbond(int _M1, int _M2, int _arm1,int _arm2){M1=_M1;        M2=_M2;   arm1=_arm1;arm2=_arm2;}
      hbond(){M1=0;      M2=0;  	arn1=0;arm2=0;}

};
#endif