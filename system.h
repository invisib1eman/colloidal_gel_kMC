//NANOROD: system.h System Class (Revision Date: Oct 27,2023)
#ifndef _SYSTEM_H
#define _SYSTEM_H
#include "header.h"
#include "xyz.h"
#include "quarternion.h"
#include "molecule.h"
#include "utils.h"
#include "grid.h"
class System
{
  public:
    vector<Molecule> M; //List of molecules
    vector<Grid> G;
    int NGRID,NGRID3;//NGRID=#cell in one direction, NGRID3=total # of cells of grid, NCELL=#particles per cell
    double NCELL;
    double GRIDL;//Grid length in one direction
    const gsl_rng_type * gsl_T;
    gsl_rng * gsl_r;
    string Description;
    int NMOL; //Number of molecules
    double arm_L=1.1;
    double bond_length=0.28;
    double bond_extension=0.1;
    double maxl2_bond=pow(bond_length+bond_extension,2);
    double cm_L=arm_L*2+bond_length;
    double max2_cm=pow(cm_L+bond_extension,2);
    double L; //Length of box
    int GSL_SEED; //Seed of random number generator
    int nsweep; //Number of MC sweeps
    double deltat; //Timestep
    double MCstep; //Step size of translation

    void ReadInput(int argc, char *argv[])
    {
        double total_time;
        
        options_description desc("Usage:\nNANOROD <options>");
    
        desc.add_options()
        ("help,h", "print usage message")
        ("NGRID,G",value<int>(&NGRID)->default_value(10),"grids(default 5)")
        ("NMOL,N", value<int>(&NMOL)->default_value(100), "#molecules (default 100)")
        ("box_length,L", value<double>(&L)->default_value(10.0), "length of box (default 10.0)")
        ("time,s", value<double>(&total_time)->default_value(100.0), "total time in tau_0 units (default 100.0)")
        ("MCstep,m", value<double>(&MCstep)->default_value(0.05), "MC step size (default 0.05)")
        ("GSL_SEED,g", value<int>(&GSL_SEED)->default_value(10), "seed for the RNG (default 10)")
        ("Description,D", value<string>(&Description)->default_value("nanorod"), "Description (default nanorod)");
;
        
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
        
        if (vm.count("help"))
        {
            cout << desc << "\n";
            exit(1);
        }
        
        gsl_rng_env_setup();
          
        gsl_T = gsl_rng_default;
        gsl_r = gsl_rng_alloc(gsl_T);
        gsl_rng_set(gsl_r, GSL_SEED);
        
        deltat=1.0/12.0*MCstep*MCstep;
        nsweep=int(ceil(total_time/deltat));
        NGRID3=NGRID*NGRID*NGRID;
        GRIDL=L/NGRID;
        NCELL=NMOL/NGRID3;
    }
    
    void Create();
    void WriteMol2(int timestep);
    void WriteDump(int timestep);
};
#endif
