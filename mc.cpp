//NANOROD: mc.cpp MC Class Function Definitions (Revision Date: October 27,2023)
#include "mc.h" 
void MC::Sweep()
{
    Mnew=S.M;
    
    double accept=0.0;
   
    int nsample=S.nsweep/100;
    if(nsample<1)
      nsample=1;
    
    energy=TotalEnergy();
    cout<<"Initial Energy="<<energy<<endl;
    
    WriteTemplate();
	LogProfile(0,accept);
    
    for(int i=1; i<=S.nsweep; i++)
    {
        time+=S.deltat;
        accept+=MoveMolecule();
        
        if(i%nsample==0)
        {
            accept/=double(nsample);
            LogProfile(i,accept);
            S.WriteMol2(i);
            S.WriteDump(i);
            accept=0.0;
        }
    }
	
   cout<<"Final Energy="<<energy<<endl;
   cout<<"Final Energy by recalculation="<<TotalEnergy()<<endl;
   S.WriteMol2(S.nsweep);
   S.WriteDump(S.nsweep);
}

void MC::WriteTemplate()
{
    ofstream out;
    out.open("_MC.log");
    out<<setw(12)<<"sweep"<<"\t"<<setw(12)<<"time"<<"\t"<<setw(12)<<"energy"<<"\t"<<setw(8)<<"accept"<<endl;
    out.close();
}

void MC::LogProfile(int i, double accept)
{
    ofstream out;
    out.open("_MC.log", ios::app);
    out<<setw(12)<<i<<"\t"<<setw(12)<<time<<"\t"<<setw(12)<<energy<<"\t"<<setw(8)<<accept<<endl;
    out.close();
}


//Returns acceptance fraction
double MC::MoveMolecule()
{
    XYZ newpos;
    double rnew, rold, delta=0.0, accept=0.0;
    double dalpha, dbeta, dgamma;
    XYZ dxyz;
    
    for(int k=0; k<S.NMOL; k++)
    {
        int i=gsl_rng_uniform_int(S.gsl_r,S.NMOL);
        
        Molecule newmolecule=S.M[i];
        newmolecule.centre=RandomTranslate(S.M[i].centre,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        newmolecule.orientation=RandomRotate(S.M[i].orientation,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        newmolecule.UpdateVertices();
        S.M[i]=newmolecule;
    }
    return accept/double(S.NMOL);
}

bool MC::Glauber(double delta, double rand)
{
    if(1.0/(exp(delta)+1.0)>rand)
      return true;
    else
      return false;
}

double MC::TotalEnergy()
{
    double totalenergy=0.0;
    double r2;
    return totalenergy;
}
