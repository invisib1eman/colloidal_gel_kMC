//NANOROD: hbond Class to manage the bond information  (Revision Date: Oct 27, 2023)
#ifndef _bond_H
#define _bond_H
#include "header.h"
class bond
{
  public:
      int P1,P2;
      void set(int _P1, int _P2){P1=P1; P2=P2;}
      bond(int _P1, int _P2){P1=_P1; P2=_P2;}
      bond(){P1=0; P2=0;}

};
#endif