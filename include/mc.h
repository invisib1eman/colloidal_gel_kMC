//NANOROD: mc.h MC Class (Revision Date: Oct 27, 2023)
#ifndef _MC_H
#define _MC_H
#include "system.h"
#include "utils.h"
#include "particle.h"
// #include "hbond.h"
#include "energies.h"
#include "grid.h"
class MC
{
    public:
        System& S;
        Energy E;
        double energy = 0.0;
        double time = 0.0;
        int N_frame;
        int nbr_g=27;//number of number grids
        //new vector for molecule
        vector<Particle> Pnew;
        // MC constructor
        MC(System& sys) : S(sys) {
            //Initialize the energy class
            E.debye_length = sys.debye_length;
            E.bjerrum_length = sys.bjerrum_length;
            E.charge = sys.charge;
            E.BoxLength = sys.BoxLength;
            E.R_hardcore = sys.R_hardcore;
            E.well_width = sys.well_width;
            E.well_edge = sys.R_hardcore + sys.well_width;
            E.well_depth = sys.well_depth;
            // Setup the R_hardcore_DH the same as the well_edge (hardcore radius)
            E.R_hardcore_DH = sys.R;
            E.cutoff_distance = sys.cutoff_distance;
            // Pass the Nframe to the MC class
            N_frame = sys.N_frame;
            E.morse_a = sys.morse_a;
            E.morse_r0 = sys.morse_r0;
            E.morse_well_depth = sys.morse_well_depth;
        }
        void WriteTemplate();
        void LogProfile(int, double );
        void Sweep();
        double MoveParticle_Single_Particle();
        double MoveParticle_Cluster_Rigid();
        double MoveParticle_Cluster_Free_Roll();
        double MoveParticle_Cluster_Alpha();
        double MoveParticle_Cluster_Alpha_Morse();
        bool Glauber(double, double);
        void Cluster_Particles();
        void Find_Neighbors();
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

