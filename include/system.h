//NANOROD: system.h System Class (Revision Date: Oct 27,2023)
#ifndef _SYSTEM_H
#define _SYSTEM_H
#include "header.h"
#include "xyz.h"
// #include "quarternion.h"
#include "utils.h"
#include "grid.h"
#include "aggregate.h"
#include "particle.h"
class System
{
  public:
    vector<Particle> P; //List of molecules
    vector<Aggregate> Ag;
    vector<Grid> G;//grid list
    int NGRID,NGRID3;//NGRID=#cell in one direction, NGRID3=total # of cells of grid, NCELL=#particles per cell
    double NCELL;
    double GRIDL;//Grid length in one direction
    const gsl_rng_type * gsl_T;
    gsl_rng * gsl_r;
    string Description;
    int NMOL; //Number of molecules
    //volume fraction 6% (core+shell)
    double R=1;//radius = 8nm
    double R_hardcore = 1.20;//hardcore radius = 1.20*8 nm
    double charge = 12.6;//charge of the molecule
    double debye_length = 0.91875;//debye length 7.35 nm at 0.9 mM salt
    double bjerrum_length = 0.1875;//bjerrum length
    double cutoff_distance = 4 * R_hardcore;//cutoff distance = 4*R_hardcore = 4*9.6 nm = 38.4 nm
    double arm_L=0.25;//arm length = 2nm
    double cm_L=(arm_L+R)*2;//the cm distance when it is possible to form a bond
    double search2_cm=pow(cm_L,2);//the cm distance when it is possible to form a bond
    double well_width = 0.05;//well width
    double L; //Length of box
    int GSL_SEED; //Seed of random number generator
    int nsweep; //Number of MC sweeps
    double deltat; //Timestep
    double MCstep; //Step size of translation
    double free_bond_freeenergy=-1;//free bond entropy
    bool freeroll = 1;//1 is true 0 is false
    bool fake_acceleration = 0;//1 is true 0 is false
    bool read_restart = 0;//1 is true 0 is false
    double tau_crawl = 100; // crawl time scale is 100*dt
    double p_crawl = (1.0/tau_crawl)/(1.0/tau_crawl+1.0); // probability of crawl
    string dump_file_name;
    string read_restart_file_name;
    string restart_file_name;
    string data_file_name;
    string log_file_name;
    void ReadInput(int argc, char *argv[])
    {
        double total_time;
        
        options_description desc("Usage:\nNANOROD <options>");
    
        desc.add_options()
        ("help,h", "print usage message")
        ("NGRID,G",value<int>(&NGRID)->default_value(25),"grids(default 10)")
        ("NMOL,N", value<int>(&NMOL)->default_value(1000), "#molecules (default 400)")
        ("box_length,L", value<double>(&L)->default_value(50.0), "length of box (default 20.0)")
        ("time,s", value<double>(&total_time)->default_value(100.0), "total time in tau_0 units (default 100.0)")
        ("MCstep,m", value<double>(&MCstep)->default_value(0.1), "MC step size (default 0.1)")// fluctuation is 0.8 nm
        ("GSL_SEED,g", value<int>(&GSL_SEED)->default_value(10), "seed for the RNG (default 10)")
        ("debye_length,d", value<double>(&debye_length)->default_value(0.91875), "debye length (default 0.91875)")
        ("Description,D", value<string>(&Description)->default_value("nanorod"), "Description (default nanorod)")
        ("well_width,w", value<double>(&well_width)->default_value(0.1), "well width (default 0.1)")
        ("freeroll,f", value<bool>(&freeroll)->default_value(1), "freeroll (default 1)")
        ("fake_acceleration,a", value<bool>(&fake_acceleration)->default_value(0), "fake acceleration (default 0)")
        ("read_restart,r", value<bool>(&read_restart)->default_value(0), "read restart (default 0)")

        ("tau_crawl,tcrawl", value<double>(&tau_crawl)->default_value(100), "crawl time scale (default 100)")
        ("dump_file_name,dump", value<string>(&dump_file_name)->default_value("dump.xyz"), "dump file name (default dump.xyz)")
        ("restart_file_name,restart", value<string>(&restart_file_name)->default_value("restart.xyz"), "restart file name (default restart.xyz)")
        ("read_restart_file_name,read_restart", value<string>(&read_restart_file_name)->default_value("restart.xyz"), "read restart file name (default restart.xyz)")
        ("data_file_name,data", value<string>(&data_file_name)->default_value("data.xyz"), "data file name (default data.xyz)")
        ("log_file_name,log", value<string>(&log_file_name)->default_value("log.txt"), "log file name (default log.txt)");
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
        // D = kT / 6πηa = 4.11×10−21/(6*pi*8*10**-9*0.8*10**-3)=3.40691*10^-11m^2/s
        //real time scale is MCstep*MCstep*(1/12)*((8*10**-9)**2)/(3.40691*10^-11)s=1.56545*10^-9s
        nsweep=int(ceil(total_time/deltat));
        NGRID3=NGRID*NGRID*NGRID;
        p_crawl = (1.0/tau_crawl)/(1.0/tau_crawl+1.0);
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
        if(GRIDL<R_hardcore*2)
        {
            cout<<"Error: Grid size too small"<<endl;
            exit(1);
        }
        NCELL=NMOL/NGRID3;
    }
    
    void Create();
    // void WriteMol2(int timestep);
    void WriteDump(int timestep);
    // void WriteBond(int timestep);
    // void WriteOrientation(int timestep);
    // void WriteGrid(int timestep);
    void UpdateGrid();
    void WriteData(int timestep);
    void CreateDump();
    void writerestart();
    void readrestart();
};
#endif
