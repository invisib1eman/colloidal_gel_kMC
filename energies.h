//MIXING: enegies.h Energy Class (Revision Date: Feb 2, 2012)
#ifndef _ENERGY_H
#define _ENERGY_H

class Energy
{
public:
    double Kd=1.0;
    double Kt=1.0;
    double Kp=1.0;
    double d0=0.3;
    double t0=0.0;
    double p0=0.0;
    double hbond(d,t,p)
    {
        return 0.5*(d-d0)*(d-d0)+0.5*(t-t0)*(t-t0)+0.5*(p-p0)*(p-p0);
    }
};
#endif

