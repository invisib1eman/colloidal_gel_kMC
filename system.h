//NANOROD: system.h System Class (Revision Date: Oct 27,2023)
#ifndef _SYSTEM_H
#define _SYSTEM_H
#include "header.h"
#include "xyz.h"
#include "utils.h"
#include "grid.h"
#include "particle.h"
class System
{
  public:
    vector<Particle> P; //List of molecules
    vector<Grid> G;//grid list
    
    int NGRID,NGRID3;//NGRID=#cell in one direction, NGRID3=total # of cells of grid, NCELL=#particles per cell
    double NCELL;
    double GRIDL;//Grid length in one direction
    const gsl_rng_type * gsl_T;
    gsl_rng * gsl_r;
    string Description;
    int NMOL; //Number of colloidal particles
    double radius=4.2;//radius is 4.2
    double L; //Length of box
    int GSL_SEED; //Seed of random number generator
    int nsweep; //Number of MC sweeps
    double deltat; //Timestep
    double MCstep; //Step size of translation
    double V_0=5;//potential well
    void ReadInput(int argc, char *argv[])
    {
        double total_time;
        
        options_description desc("Usage:\nNANOROD <options>");
    
        desc.add_options()
        ("help,h", "print usage message")
        ("NGRID,G",value<int>(&NGRID)->default_value(10),"grids(default 10)")
        ("NMOL,N", value<int>(&NMOL)->default_value(100), "#molecules (default 400)")
        ("box_length,L", value<double>(&L)->default_value(677), "length of box (default 20.0)")
        ("time,s", value<double>(&total_time)->default_value(100.0), "total time in tau_0 units (default 100.0)")
        ("MCstep,m", value<double>(&MCstep)->default_value(0.1), "MC step size (default 0.1)")
        ("GSL_SEED,g", value<int>(&GSL_SEED)->default_value(10), "seed for the RNG (default 10)")
        ("Potential_well,v", value<double>(&V_0)->default_value(10), "Potential well for Morse (default 10)")

        ("Description,D", value<string>(&Description)->default_value("nanorod"), "Description (default nanorod)");
        cout<<NMOL<<endl;

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
        
    try
    {
      G.reserve(NGRID3);
    }
    catch (int e)
    {
	cout<<"Memory issues in cell list allocation.. exiting"<<endl;
	exit(1);
    }
        GRIDL=L/NGRID;
        if(GRIDL<radius*2)
        {
            cout<<"Error: Grid size too small"<<endl;
            exit(1);
        }
        NCELL=NMOL/NGRID3;
    }
    
    void Create();
    void WriteDump(int timestep);
    void WriteGrid(int timestep);
    void UpdateGrid();
};
#endif
