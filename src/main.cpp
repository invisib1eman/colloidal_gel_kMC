//NANOROD: main.cpp Main Program (Revision Date: Oct 27, 2023)
#include "header.h"
#include "system.h"
#include "mc.h"

int main(int argc, char *argv[]) {
    time_t start, end;
    time(&start);
    System sys;
    sys.ReadInput(argc,argv);
    if(sys.read_restart==1) {
        sys.readrestart();
    } else {
        sys.Create();
    }
    MC mc(sys);
    mc.N_frame = sys.N_frame;
    // Pass parameters from sys to energy class
    mc.E.debye_length = sys.debye_length;
    mc.E.bjerrum_length = sys.bjerrum_length;
    mc.E.charge = sys.charge;
    mc.E.L = sys.L;
    mc.E.R_hardcore = sys.R_hardcore;
    mc.E.well_width = sys.well_width;
    mc.E.well_edge = sys.R_hardcore + sys.well_width;
    mc.E.R_hardcore_DH = mc.E.well_edge;
    mc.Sweep();
    sys.writerestart();
    //Finalize Random number
    gsl_rng_free(sys.gsl_r);
    time(&end);
    cout<<"Time Elapsed="<<difftime(end,start)<<endl;
    return 0;
}
