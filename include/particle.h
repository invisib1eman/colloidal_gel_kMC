#ifndef _PARTICLE_H
#define _PARTICLE_H
#include "header.h"
#include "xyz.h"
//#include "quaternion.h"
#include "bond.h"
class Particle
{
public:
    int P_ID;//particle ID
    char Type; //type of particle
    int gID;//grid id
    XYZ position; //Coordinates of vertices
    int nbonds;//number of bonds, 0 or 1 (only for ver)
    int A_ID;//aggregate id
    set<int> nbr_particles;
    
    
        
    
    
    
    Particle() //HARD CODE BASED ON VERTEX INFORMATION,unit length=1nm,arm length=1.1nm
    {
        
        Type='O';
        gID=0;
        position.set(0.0,0.0,0.0);
        P_ID=0;
        A_ID=0;
    }
};
#endif