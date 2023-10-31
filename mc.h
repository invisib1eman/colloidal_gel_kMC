//NANOROD: mc.h MC Class (Revision Date: Oct 27, 2023)
#ifndef _MC_H
#define _MC_H
#include "system.h"
#include "utils.h"
#include "molecule.h"
#include "energies.h"
class MC
{
    public:
        System S;
        double energy, time;
        //new vector for molecule
        vector<Molecule> Mnew;
        MC(){energy=0.0; time=0.0; }
        void WriteTemplate();
        void LogProfile(int, double );
        void Sweep();
        double MoveMolecule();
        bool Glauber(double, double);
        double TotalEnergy();
};
#endif

