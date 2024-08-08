//NANOROD: mc.h MC Class (Revision Date: Oct 27, 2023)
#ifndef _MC_H
#define _MC_H
#include "system.h"
#include "utils.h"
#include "particle.h"
#include "energies.h"
#include "grid.h"
class MC
{
    public:
        System S;
        Energy E;
        double energy, time;
        int nbr_g=27;//number of number grids
        MC(){energy=0.0; time=0.0; }
        void WriteTemplate();
        void LogProfile(int, double );
        void Sweep();
        double MoveParticle();
        bool Glauber(double, double);
        double Energy();
        
        double WriteEnergy(int timestep);
};
#endif

