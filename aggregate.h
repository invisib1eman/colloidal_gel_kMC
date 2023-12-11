//KMC: aggregate.h Class (Added: Dec 24, 2010)
//MIXING: enegies.h Energy Class (Revision Date: Feb 2, 2012)
#ifndef _AGGREGATE_H
#define _AGGREGATE_H
#include "utils.h"
#include "molecule.h"
#include "hbond.h"
#include "xyz.h"
#include "quarternion.h"
class Aggregate
{
public:
    int n; //number of particles
    list<int> plist; //List of particles' molecular id
    double rg2; //radius of gyration squared;
    Aggregate(){n=0; rg2=0.0;}
};
#endif