//KMC: aggregate.h Class (Added: Dec 24, 2010)
//MIXING: enegies.h Energy Class (Revision Date: Feb 2, 2012)
#ifndef _AGGREGATE_H
#define _AGGREGATE_H
#include "utils.h"
#include "bond.h"
#include "particle.h"
#include "xyz.h"
#include "quaternion.h"
class Aggregate
{
public:
    int n; //number of particles
    vector<int> plist; //List of particles' particle id
    double Rg; //radius of gyration squared;
    quaternion Orientation; //orientation of the aggregate
    
    XYZ cm; //center of mass
    Aggregate(){n=0; Rg=0.0; cm.x=0.0; cm.y=0.0; cm.z=0.0; plist.clear();}
    void UpdateOrientation();
};
#endif