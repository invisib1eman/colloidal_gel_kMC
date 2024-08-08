#ifndef _PARTICLE_H
#define _PARTICLE_H
#include "header.h"
#include "xyz.h"
class Particle
{
public:
    int P_ID;//particle ID
    char Type; //type of particle
    int gID;//grid id
    XYZ position; //Coordinates of vertices
    double radius;//radius of the particle
    
        
    
    
    
    Particle() //HARD CODE BASED ON VERTEX INFORMATION,unit length=1nm,arm length=1.1nm
    {
        
        Type='O';
        gID=0;
        position.set(0.0,0.0,0.0);
        radius=4.2;
        
    }
};
#endif