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
    sys.WriteMol2(0);
	
    MC mc;
    mc.S=sys;
    mc.Sweep();
    
	//Finalize Random number
	gsl_rng_free(sys.gsl_r);
	time(&end);
	cout<<"Time Elapsed="<<difftime(end,start)<<endl;
    return 0;
}
