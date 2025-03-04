//KMC: grid.h Grid Class (Revision Date: Dec 22, 2010)
#ifndef _GRID_H
#define _GRID_H
#include "header.h"
#include "xyz.h"
class Grid
{
public:
    XYZ cm;
    int g_index;
    int n;//Number of particles
    list<int> plist; 
    vector<int> nbr; //List of neighboring cells including itself: NOT IMPLEMENTED IN GRID BUT IN CELLLIST
    Grid(){cm.x=0.0; cm.y=0.0; cm.z=0.0; n=0;}
};
#endif
