#ifndef _PARTICLE_H
#define _PARTICLE_H
#include "header.h"
#include "xyz.h"
//#include "quarternion.h"
#include "hbond.h"
class Particle
{
public:
    int P_ID;//particle ID
    char Type; //type of particle, can be O(center) or A B C D or E F G H(corresponding to bonded A B C D)
    int gID;//grid id
    XYZ position; //Coordinates of vertices
    int nbonds;//number of bonds, 0 or 1 (only for ver)
    
    
        
    
    vector<hbond> hbond_list;//list of hbonded neighbors and information, size be 0 or 1 (only for ver)
    
    
    Particle() //HARD CODE BASED ON VERTEX INFORMATION,unit length=1nm,arm length=1.1nm
    {
        
        nbonds=0;
        Type='O';
        gID=0;
        position.set(0.0,0.0,0.0);
        hbond_list.clear();
        
    }
};
#endif