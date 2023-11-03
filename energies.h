//MIXING: enegies.h Energy Class (Revision Date: Feb 2, 2012)
#ifndef _ENERGY_H
#define _ENERGY_H
#include "system.h"
#include "utils.h"
#include "molecule.h"
#include "hbond.h"
#include "xyz.h"
#include "quarternion.h"
class Energy
{
public:
    double Kr=10.0;
    double Kalpha=100.0;
    double Kxhi=100.0;
    double r0=0.28;
    double alpha0=0.0;
    double beta0=0;
    double xhi0=0.0;
    double L=10;//box size
    double hbonde(Molecule M1,Molecule M2,int n)
    {
        int armindex1=M1.hbond_list[n].arm1;
        XYZ arm1=M1.ver[armindex1];
        XYZ centre1=M1.centre;
        XYZ neighborarm1=M1.ver[armindex1+1-2*(armindex1%2)];
        int armindex2=M1.hbond_list[n].arm2;
        XYZ arm2=M2.ver[armindex2];
        XYZ centre2=M2.centre;
        XYZ neighborarm2=M2.ver[armindex2+1-2*(armindex2%2)];
        double r=sqrt(min_d2(arm1,arm2,L));
        double alpha=angle_vectors((arm1-centre1),(arm2-arm1));
        double beta=angle_vectors((arm2-centre2),(arm1-arm2));
        double xhi=dihedral_vectors((centre1-neighborarm1),(centre2-centre1),(neighborarm2-centre2));
        return 0.5*Kr*(r-r0)*(r-r0)+0.5*Kalpha*(alpha-alpha0)*(alpha-alpha0)+0.5*Kalpha*(beta-beta0)*(beta-beta0)+0.5*Kxhi*(xhi-xhi0)*(xhi-xhi0);
    }
};
#endif

