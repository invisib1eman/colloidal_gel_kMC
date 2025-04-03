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
    // mode: "cluster_rigid" or "cluster_free_roll" or "single_particle" or "cluster_alpha" or "cluster_alpha_morse"
    string mode;
    double alpha=0.33;
    double morse_range=1.0;
    double morse_a=0.5;
    double morse_r0=2.0;
    double morse_well_depth=-5;
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
    // volume fraction (core+shell) 
    // define contants (not change as parameters)
    double R = 1;//radius = 8nm, the semiconductor core radius and unit length scale
    double R_hardcore = 1.20;//hardcore radius = 1.20*8 nm
    double charge = 12.6;//charge of the molecule
    double debye_length = 0.91875;//debye length 7.35 nm at 0.9 mM salt
    double bjerrum_length = 0.1875;//bjerrum length
    double well_width = 0.05;//well width
    double well_edge = R_hardcore + well_width; //well edge = 10nm
    double well_depth = -10000;
    double cm_L = well_edge * 2;//the cm distance when it is possible to form a bond
    double search2_cm = pow(cm_L,2);//the cm distance when it is possible to form a bond
    double cutoff_distance; //cutoff distance = 2 * well_edge + 4 * debye_length
    double BoxLength; //Length of box
    int GSL_SEED; //Seed of random number generator
    int nsweep; //Number of MC sweeps
    double deltat; //Timestep
    double MCstep; //Step size of translation
    int N_frame; //Number of frames
    double free_bond_freeenergy=-1;//free bond entropy
    bool fake_acceleration = 0;//1 is true 0 is false
    bool read_restart = 0;//1 is true 0 is false
    bool read_xyz = 0;//1 is true 0 is false
    double tau_crawl = 100; // crawl time scale is 100*dt
    double p_crawl = (1.0/tau_crawl)/(1.0/tau_crawl+1.0); // probability of crawl
    string dump_file_name;
    string read_restart_file_name;
    string restart_file_name;
    string data_file_name;
    string log_file_name;
    string read_xyz_file_name;
    void ReadInput(int argc, char *argv[])
    {
        double total_time;
        
        options_description desc("Usage:\nNANOROD <options>");
    
        desc.add_options()
        ("help,h", "print usage message")
        ("NGRID,G",value<int>(&NGRID)->default_value(20),"grids(default 20)")
        ("NMOL,N", value<int>(&NMOL)->default_value(1000), "#molecules (default 1000)")
        ("box_length,L", value<double>(&BoxLength)->default_value(50.0), "length of box (default 50.0)")
        ("time,s", value<double>(&total_time)->default_value(100.0), "total time in tau_0 units (default 100.0)")
        ("MCstep,m", value<double>(&MCstep)->default_value(0.1), "MC step size (default 0.1)")// fluctuation is 0.8 nm
        ("GSL_SEED,g", value<int>(&GSL_SEED)->default_value(10), "seed for the RNG (default 10)")
        ("debye_length,d", value<double>(&debye_length)->default_value(0.91875), "debye length (default 0.91875)")
        ("Description,D", value<string>(&Description)->default_value("colloidgel"), "Description (default colloidgel)")
        ("well_width,w", value<double>(&well_width)->default_value(0.05), "well width (default 0.05)")
        ("well_depth,W", value<double>(&well_depth)->default_value(-10000), "well depth (default -10000)")
        ("mode", value<string>(&mode)->default_value("cluster_free_roll"), "mode (default cluster_free_roll) (choices: cluster_free_roll, cluster_rigid, single_particle, cluster_alpha, cluster_alpha_morse)")
        ("alpha", value<double>(&alpha)->default_value(0.0), "alpha (default 0.33) (D~N^(-alpha))")
        ("morse_range,M", value<double>(&morse_range)->default_value(1.0), "morse range (default 1.0)")
        ("morse_a,A", value<double>(&morse_a)->default_value(0.5), "morse a (default 0.5)")
        ("morse_r0,R", value<double>(&morse_r0)->default_value(1.0), "morse r0 (default 1.0)")
        ("morse_well_depth,W", value<double>(&morse_well_depth)->default_value(-5), "morse well depth (default -5)")
        ("fake_acceleration,a", value<bool>(&fake_acceleration)->default_value(0), "fake acceleration (default 0)")
        ("read_restart,r", value<bool>(&read_restart)->default_value(0), "read restart (default 0)")
        ("tcrawl", value<double>(&tau_crawl)->default_value(100), "crawl time scale (default 100)")
        ("dump", value<string>(&dump_file_name)->default_value("a.lammpsdump"), "dump file name (default a.lammpsdump)")
        ("restart", value<string>(&restart_file_name)->default_value("a.restart"), "restart file name (default a.restart)")
        ("read_restart_file", value<string>(&read_restart_file_name)->default_value("a.restart"), "read restart file name (default a.restart)")
        ("data", value<string>(&data_file_name)->default_value("a.data"), "data file name (default a.data)")
        ("log,l", value<string>(&log_file_name)->default_value("a.log"), "log file name (default a.log)")
        ("read_xyz,x", value<bool>(&read_xyz)->default_value(0), "read xyz (default 0)")
        ("read_xyz_file", value<string>(&read_xyz_file_name)->default_value("a.xyz"), "read xyz file name (default a.xyz)")
        ("nframe,n", value<int>(&N_frame)->default_value(1000), "Number of frames (default 1000)");
        // define the input machine
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
        
        if (vm.count("help"))
        {
            cout << desc << "\n";
            exit(1);
        }
        // define the random number generator
        gsl_rng_env_setup();
        gsl_T = gsl_rng_default;
        gsl_r = gsl_rng_alloc(gsl_T);
        gsl_rng_set(gsl_r, GSL_SEED);
        // Calculate variables from input
        // calculate the well edge
        well_edge = R_hardcore + well_width;
        // calculate the cm distance
        cm_L = well_edge * 2;
        // calculate the search2_cm
        search2_cm = pow(cm_L,2);
        // calculate the cutoff distance
        cutoff_distance = 2 * well_edge + 3 * debye_length;
        // calculate the timestep
        deltat=1.0/12.0*MCstep*MCstep*(1-p_crawl);
        morse_a = -log(1-sqrt(1-exp(-morse_range/morse_a)))/morse_range;
        // D = kT / 6πηa = 4.11×10−21/(6*pi*8*10**-9*0.8*10**-3)=3.40691*10^-11m^2/s in the mixed solvent DMF/Ethylene glycol with viscosity 0.8 mPa*s.
        // real time scale of each MC step is MCstep*MCstep*(1/12)*((8*10**-9)**2)/(3.40691*10^-11)s=1.56545*10^-9s
        // calculate the number of sweeps
        if (mode == "single_particle" or mode == "cluster_rigid")
        {
            nsweep=int(ceil(total_time/deltat));
        }
        if (mode == "cluster_free_roll" or mode == "cluster_alpha" or mode == "cluster_alpha_morse")
        {
            p_crawl = (1.0/tau_crawl)/(1.0/tau_crawl+1.0);
            nsweep=int(ceil(total_time/deltat));
        }
        NGRID3=NGRID*NGRID*NGRID;
        //reserve memory for grid list
        try
        {
        G.reserve(NGRID3);
        }
        catch (int e)
        {
        cout<<"Memory issues in cell list allocation.. exiting"<<endl;
        exit(1);
        }
        GRIDL=BoxLength/NGRID;
        if(GRIDL < cutoff_distance)
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
    void WriteRestart();
    void ReadRestart();
    void ReadXYZ();
};
#endif
