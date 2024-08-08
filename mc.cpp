//NANOROD: mc.cpp MC Class Function Definitions (Revision Date: October 27,2023)
#include "mc.h" 
void MC::Sweep()
{
    E.L=S.L;
    S.WriteDump(0);
    double accept=0.0;
   
    int nsample=S.nsweep/100;
    if(nsample<1)
      nsample=1;
    
    energy=Energy();
   
    
    WriteTemplate();
	LogProfile(0,accept);
    //Preparation
    int N_preparation=S.nsweep/100;
    E.V_0=0;
    
    for(int i=1; i<=N_preparation; i++)
    {
        time+=S.deltat;
        accept+=MoveParticle();
        
        if(i%nsample==0)
        {
            
            LogProfile(i,accept);
            WriteEnergy(i);
            S.WriteDump(i);
            S.WriteGrid(i);
            accept=0.0;
        }
    }
    //Turn on attraction
    E.V_0=S.V_0;
	for(int i=N_preparation; i<=S.nsweep; i++)
    {
        time+=S.deltat;
        accept+=MoveParticle();
        
        if(i%nsample==0)
        {
            
            LogProfile(i,accept);
            WriteEnergy(i);
            S.WriteDump(i);
            S.WriteGrid(i);
            accept=0.0;
        }
    }
   
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
double MC::MoveParticle()
{
    //S.UpdateGrid();
    //check if the total particle in the celllists add to N
    int total_particles=0;

    for(int i=0;i<S.NGRID3;i++)
    {
        /*if(S.G[i].nbr.size()!=nbr_g)
        {
            cout<<"Neighbor Number wrong"<<endl;
            exit(0);
        }*/
        if(S.G[i].n!=S.G[i].plist.size())
        {
            cout<<"Number wrong"<<endl;
            exit(0);
        }
        total_particles+=S.G[i].n;
    }
    if(total_particles!=S.NMOL)
    {
            cout<<"Total number wrong"<<endl;
            exit(0);
    }
    ofstream out;
    out.open("events.txt",ios::app);
    ofstream out2;
    out2.open("Error.txt",ios::app);
    //ofstream out2;
    //out2.open("bond_energy.txt",ios::app);
    double accept=0.0;//accept events
    
    int index;
    double delta;
    
    
    //trial move
    for(int i=0; i<S.NMOL; i++)
    {
        index=gsl_rng_uniform_int(S.gsl_r,S.NMOL);
        
        delta=0;//energy difference
        //int index=sequence_M[i];//index of the current trial particle
        
        
        
       
        Particle newparticle=S.P[index];
        newparticle.position=RandomTranslate(S.P[index].position,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        XYZ image_center=image(newparticle.position,S.L);
        int new_gID=GridIndex_xyz(image_center,S.NGRID,S.GRIDL,S.L);
        newparticle.gID=new_gID;
        
        
        //add wca energy
        
        //add new wca
        for(int l=0;l<nbr_g;l++)
        {
            //ID of new neighboring grids
            int temp_gID=S.G[new_gID].nbr[l];
            if(temp_gID>=S.NGRID3)
            {
                cout<<new_gID<<"Error: wrong grid id"<<endl;
                cout<<temp_gID<<endl;
                exit(0);
            }
            if(S.G[temp_gID].plist.empty())
            {
                continue;
            }
            //iterate about particles in the neighborlist
            list<int>::iterator it;
            for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
            {
                int j=*it;
                if(j!=index)
                {
                    delta+=(E.WCA(newparticle,S.P[j])+E.Morse(newparticle,S.P[j]));
                    
                }
            }
        }
        //subtract old wca
        int old_gID=S.P[index].gID;
        for(int l=0;l<nbr_g;l++)
        {
            //ID of new neighboring grids
            int temp_gID=S.G[old_gID].nbr[l];
            list<int>::iterator it;
            for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
            {
                int j=*it;
                if(j!=index)
                {
                    
                    delta-=(E.WCA(S.P[index],S.P[j])+E.Morse(S.P[index],S.P[j]));


                }
            }
        }
        if(Glauber(delta,gsl_rng_uniform(S.gsl_r)))
        {
            //Accept move
            S.P[index]=newparticle;
    
            accept+=1.0;
            energy+=delta;

            //update grid list
            if(new_gID!=old_gID)
            {
                S.G[old_gID].n-=1;
                S.G[old_gID].plist.remove(S.P[index].P_ID);
                
                S.G[new_gID].n+=1;
                S.G[new_gID].plist.push_back(S.P[index].P_ID);
            }
        
            
        }
        
        
    }
    
    return accept/double(S.NMOL);
}

bool MC::Glauber(double delta, double rand)
{
    if (delta>100000)
        return false;
    else
    {
        if(1.0/(exp(delta)+1.0)>rand)
            return true;
        else
            return false;
    }
    
}

double MC::Energy()
{
    double ewca=0;
    list<int>::iterator it;
    for(int i=0;i<S.NMOL;i++)
    {
        int old_gID=S.P[i].gID;
        for(int l=0;l<nbr_g;l++)
        {
            //ID of new neighboring grids
            int temp_gID=S.G[old_gID].nbr[l];
            
            for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
            {
                int j=*it;
                if(j>i)//avoid double counting
                {
                    //image distance
                    
                    ewca+=E.WCA(S.P[i],S.P[j])+E.Morse(S.P[i],S.P[j]);
                    


                }
            }
        }
    }
    return ewca;
}






double MC::WriteEnergy(int timestep)
{
    ofstream out;
    out.open("energy.txt",ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<"Energy"<<endl;
    out<<Energy()<<endl;
    out.close();
    return 0;
}
