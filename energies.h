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
    double Kr=20.0;
    double Kalpha=300.0;
    double Kxhi=50.0;
    double r0=0.28;
    double fene_R0=0.58;//max_bond_extension
    double alpha0=0.0;
    double beta0=0;
    double xhi0=0.0;
    double L;//box size
    double sigma=1.8;//2*coresize
    double sigmaarm=0.28;//FENE hardcore diameter
    double epsilon=20;//LJ strength
    double feneepsilon=1;
    double hbonde(Molecule M1,Molecule M2,int n)
    {
        int armindex1=M1.hbond_list[n].arm1;
        XYZ arm1=image(M1.ver[armindex1],L);
        XYZ centre1=image(M1.centre,L);
        XYZ neighborarm1=image(M1.ver[neighborarm(armindex1)],L);
        int armindex2=M1.hbond_list[n].arm2;
        XYZ arm2=image(M2.ver[armindex2],L);
        XYZ centre2=image(M2.centre,L);
        XYZ neighborarm2=image(M2.ver[neighborarm(armindex2)],L);
        double r=sqrt(min_d2(arm1,arm2,L));
        double alpha=angle_vectors(real_vector(arm1-centre1,L),real_vector(arm2-arm1,L));
        double beta=angle_vectors(real_vector(arm2-centre2,L),real_vector(arm1-arm2,L));
        double xhi=dihedral_vectors(real_vector(neighborarm1-centre1,L),real_vector(centre2-centre1,L),real_vector(neighborarm2-centre2,L));
        return 0.5*Kr*(r-r0)*(r-r0)+0.5*Kxhi*(xhi-xhi0)*(xhi-xhi0)+0.5*Kalpha*(alpha-alpha0)*(alpha-alpha0)+0.5*Kalpha*(beta-beta0)*(beta-beta0);
        
    };
    double hbonde_fene(Molecule M1,Molecule M2,int armindex1,int armindex2)
    {
        
        XYZ arm1=image(M1.ver[armindex1],L);
        
        
        XYZ arm2=image(M2.ver[armindex2],L);
        double wca=0;
        double fene=0;
        double r=sqrt(min_d2(arm1,arm2,L));
        double r2=r*r;
        double nr2=r2/(sigmaarm*sigmaarm);
        double nr2cut=1.12*1.12;
        double nr6=nr2*nr2*nr2;
        double nr6cut=nr2cut*nr2cut*nr2cut;
        if(nr2<nr2cut)//2^(1/6)=1.12
        {
            
            wca=4*feneepsilon*(1/(nr6*nr6)-(1/nr6))-4*feneepsilon*(1/(nr6cut*nr6cut)-(1/nr6cut));
        }
        else
        {
            
            
            wca=0;
        }
        if (r>fene_R0)
        {
            fene=10000000;//not allowed to extend more than fene_R0
        }
        else
        {
            
            fene=-0.5*Kr*fene_R0*fene_R0*log(1-(r/fene_R0)*(r/fene_R0));
        }
        return wca+fene;
        
    };
    double hbonde_angle(Molecule M1,Molecule M2,int armindex1,int armindex2)
    {
        
        XYZ arm1=image(M1.ver[armindex1],L);
        XYZ centre1=image(M1.centre,L);
        
        
        XYZ arm2=image(M2.ver[armindex2],L);
        XYZ centre2=image(M2.centre,L);
        
       
        
    
        
        double alpha=angle_vectors(real_vector(arm1-centre1,L),real_vector(arm2-arm1,L));
        double beta=angle_vectors(real_vector(arm2-centre2,L),real_vector(arm1-arm2,L));
            
        return 0.5*Kalpha*(alpha-alpha0)*(alpha-alpha0)+0.5*Kalpha*(beta-beta0)*(beta-beta0);
        
    }
    double hbonde_dihedral(Molecule M1,Molecule M2,int armindex1,int armindex2)
    {
        
        XYZ arm1=image(M1.ver[armindex1],L);
        XYZ centre1=image(M1.centre,L);
        XYZ neighborarm1=image(M1.ver[neighborarm(armindex1)],L);
        
        XYZ arm2=image(M2.ver[armindex2],L);
        XYZ centre2=image(M2.centre,L);
        XYZ neighborarm2=image(M2.ver[neighborarm(armindex2)],L);
        double xhi=dihedral_vectors(real_vector(neighborarm1-centre1,L),real_vector(centre2-centre1,L),real_vector(neighborarm2-centre2,L));
        return 0.5*Kxhi*(xhi-xhi0)*(xhi-xhi0);
        
    }
    double LJ(double r2)
    {
        double nr2=r2/(sigma*sigma);
        double nr6=nr2*nr2*nr2;
        return 4*epsilon*(1/(nr6*nr6)-(1/nr6));

    };
    double WCA(double r2)
    {
        double nr2=r2/(sigma*sigma);
        double nr2cut=1.12*1.12;
        double nr6=nr2*nr2*nr2;
        double nr6cut=nr2cut*nr2cut*nr2cut;
        if(nr2<nr2cut)//2^(1/6)=1.12
        {
            
            return 4*epsilon*(1/(nr6*nr6)-(1/nr6))-4*epsilon*(1/(nr6cut*nr6cut)-(1/nr6cut));
        }
        else
        {
            
            
            return 0;
        }
    };
};
#endif

