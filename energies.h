//MIXING: enegies.h Energy Class (Revision Date: Feb 2, 2012)
#ifndef _ENERGY_H
#define _ENERGY_H
#include "system.h"
#include "utils.h"
#include "xyz.h"
#include "particle.h"
class Energy
{
public:
    double L;//box size
    double epsilon=1;//LJ strength
    double V_0;//Depth of Morse potential
    
    double LJ(Particle P1,Particle P2)
    {
        double r2=min_d2(P1.position,P2.position,L);
        double sigma=P1.radius+P2.radius;
        double nr2=r2/(sigma*sigma);
        double nr6=nr2*nr2*nr2;
        if(nr2<nr2cut)//2^(1/6)=1.12
        {
            
            return 4*epsilon*(1/(nr6*nr6)-(1/nr6));
        }
        else
        {
            
            
            return 0;
        }
    };       
    double hardcore(Particle P1,Particle P2)
    {
        double r2=min_d2(P1.position,P2.position,L);
        double sigma=P1.radius+P2.radius;
        if (r2<sigma)
            return 100001;
        else
            return 0;
            
            
           
    };  
    
    
    double WCA(Particle P1,Particle P2)
    {
        double r2=min_d2(P1.position,P2.position,L);
        double sigma=P1.radius+P2.radius;
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
    double Morse(Particle P1,Particle P2)
    {
        double r=sqrt(min_d2(P1.position,P2.position,L));
        double a1=P1.radius;
        double a2=P2.radius;
        double kappa=sqrt(3);
        return V_0*(1-exp(-kappa*(r-(a1+a2))))*(1-exp(-kappa*(r-(a1+a2))));
    };
};
#endif

