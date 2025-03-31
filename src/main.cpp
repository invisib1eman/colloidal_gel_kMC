//NANOROD: main.cpp Main Program (Revision Date: Oct 27, 2023)
#include "header.h"
#include "system.h"
#include "mc.h"

int main(int argc, char *argv[]) {
    // Start timing
    time_t start, end;
    time(&start);
    // Initialize the system, read input, and create the system
    System sys;
    sys.ReadInput(argc,argv);
    // If read_restart is true, read the restart file
    if(sys.read_restart==1) {
        sys.ReadRestart();
    } else {
        sys.Create();
    }
    if(sys.read_xyz==1) {
        sys.ReadXYZ();
    }
    // Initialize the MC class
    MC mc(sys);
    // Run the MC sweeps
    mc.Sweep();
    // Write the restart file
    sys.WriteRestart();
    //Finalize Random number
    gsl_rng_free(sys.gsl_r);
    // End timing
    time(&end);
    // Print the time elapsed
    cout<<"Time Elapsed="<<difftime(end,start)<<endl;
    // Return 0
    return 0;
}
