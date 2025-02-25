//NANOROD: mc.h MC Class (Revision Date: Oct 27, 2023)
#ifndef _MC_H
#define _MC_H
#include "system.h"
#include "utils.h"
#include "particle.h"
// #include "hbond.h"
#include "energies.h"
#include "grid.h"
#include "aggregate.h"
class MC
{
    public:
        System S;
        Energy E;
        double energy, time;
        int nbr_g=27;//number of number grids
        //new vector for molecule
        vector<Particle> Pnew;
        MC(){energy=0.0; time=0.0; }
        void WriteTemplate();
        void LogProfile(int, double );
        void Sweep();
        double MoveParticle();
        bool Glauber(double, double);
        // bool Arrhenius(double A,double delta, double rand);
        // double WCAEnergy();
        // double FENE_energy();
        // double Angle_energy();
        // double Dihedral_energy();
        // double bond_energy();
        // double bond_freeze_freenerngy();
        // double TotalEnergy();
        // double WriteEnergy(int timestep);
};
#endif

