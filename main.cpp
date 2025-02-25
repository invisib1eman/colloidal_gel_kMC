//NANOROD: main.cpp Main Program (Revision Date: Oct 27, 2023)
#include "header.h"
#include "system.h"
#include "mc.h"


int main(int argc, char *argv[])
{
    
    
	time_t start, end;
	time(&start);
	
    System sys;
    sys.ReadInput(argc,argv);
    
    sys.Create();
    // sys.WriteMol2(0);
	
    MC mc;
    mc.S=sys;
    mc.E.debye_length = sys.debye_length;
    mc.E.bjerrum_length = sys.bjerrum_length;
    mc.E.charge = sys.charge;
    mc.E.L = sys.L;
    mc.E.R_hardcore = sys.R_hardcore;
    mc.E.well_width = sys.well_width;
    mc.E.well_edge = sys.R_hardcore + sys.well_width;
    mc.Sweep();
    
	//Finalize Random number
	gsl_rng_free(sys.gsl_r);
	time(&end);
	cout<<"Time Elapsed="<<difftime(end,start)<<endl;
    return 0;
    /*
   XYZ a=quarterrotation(XYZ(1.0,0.0,-0.47),quarternion(0.899125,0.28284,0.290824,-0.164306));
   cout<<a.x<<endl;
   cout<<a.y<<endl;
   cout<<a.z<<endl;
   return 0;*/
}
