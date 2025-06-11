//NANOROD: mc.h MC Class (Revision Date: Oct 27, 2023)
#ifndef _MC_H
#define _MC_H
#include "system.h"
#include "utils.h"
#include "particle.h"
// #include "hbond.h"
#include "energies.h"
#include "grid.h"
#include "quaternion.h"
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
            E.cm_L = sys.cm_L;
            E.well_depth = sys.well_depth;
            // Setup the R_hardcore_DH the same as the R(core radius)
            E.R_hardcore_DH = sys.R;
            E.cutoff_distance = sys.cutoff_distance;
            // Pass the Nframe to the MC class
            N_frame = sys.N_frame;
            E.morse_a = sys.morse_a;
            E.morse_r0 = sys.morse_r0;
            E.morse_well_depth = sys.morse_well_depth;
            E.yukawa_e = sys.yukawa_e;
            E.yukawa_debye_length = sys.yukawa_debye_length;
            double saturation_charge = 4*E.R_hardcore_DH*(1+E.R_hardcore_DH/E.debye_length)/E.bjerrum_length;
            if (E.charge > saturation_charge)
            {
                E.charge = saturation_charge;
                cout << "Saturation charge: " << E.charge << endl;
            }
        }
        void WriteTemplate();
        void LogProfile(int, double );
        void Sweep();
        double MoveParticle_Single_Particle();
        double MoveParticle_Cluster_Rigid();
        double MoveParticle_Cluster_Free_Roll();
        double MoveParticle_Cluster_Free_Roll_Yukawa();
        double MoveParticle_Cluster_Free_Roll_Rotation();
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

