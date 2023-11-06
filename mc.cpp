//NANOROD: mc.cpp MC Class Function Definitions (Revision Date: October 27,2023)
#include "mc.h" 
void MC::Sweep()
{
    Mnew=S.M;
    S.WriteDump(0);
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
    
    double accept=0.0;//accept events
    
    //vector<int> sequence_M=generateRandom(S.NMOL);//randomly generate a sequence of the N molecules;
    for(int i=0; i<S.NMOL; i++)
    {
        int index=gsl_rng_uniform_int(S.gsl_r,S.NMOL);
        
        double delta=0;//energy difference
        //int index=sequence_M[i];//index of the current trial molecule
        Molecule newmolecule=S.M[index];
        
        newmolecule.centre=RandomTranslate(S.M[index].centre,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        

        newmolecule.orientation=RandomRotate(S.M[index].orientation,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        newmolecule.UpdateVertices();
        
        if (newmolecule.nbonds>0)   //calculate hbond_energy
        {
            
            for(int n=0; n<newmolecule.nbonds; n++)
            {
                delta+=E.hbonde(newmolecule,S.M[newmolecule.hbond_list[n].M2],n)-E.hbonde(S.M[index],S.M[newmolecule.hbond_list[n].M2],n);
                //cout<<delta<<endl;
            }
        }
        if(Glauber(delta,gsl_rng_uniform(S.gsl_r)))
        {
            //Accept move
            S.M[index]=newmolecule;
            cout<<"accept"<<endl;
            accept+=1.0;
            energy+=delta;
        //find new bonding possibility
        int new_hbond=-1;
        vector<hbond> new_hbondlist;
	    double r2_newbond=S.L*S.L;
        for (int j=0;j<S.NMOL;j++)//loop over all molecules;lately can be improved by celllist
        {   
            if (j!=index)
            {
                double r2_cm=min_d2(newmolecule.centre,S.M[j].centre,S.L);
                if (r2_cm<S.max2_cm)//check if the cm distance is in the range, can save time
                {
                    for(int k=0;k<newmolecule.N_VER;k++)
                    {
                        for(int l=0;l<S.M[j].N_VER;l++)
                        {
                            double r2_arms=min_d2(newmolecule.ver[k],S.M[j].ver[l],S.L);
                            if (r2_arms<S.maxl2_bond&&newmolecule.vertype[k]=='A'&&S.M[j].vertype[l]=='A')
                            {
                                
                                new_hbondlist.push_back(hbond(index,j,k,l));
                                new_hbond=1;
                            }
                        }
                    }
                    
                }
            }
        }
        
            //form bonds if accept
            if (new_hbond!=-1)
            {
                for(int m=0;m<new_hbondlist.size();m++) 
                {
                    S.M[index].hbond_list.push_back(new_hbondlist[m]);
                    S.M[index].vertype[new_hbondlist[m].arm1]='I';
                    S.M[index].nbonds+=1;
                    S.M[new_hbondlist[m].M2].hbond_list.push_back(hbond(new_hbondlist[m].M2,new_hbondlist[m].M1,new_hbondlist[m].arm2,new_hbondlist[m].arm1));
                    S.M[new_hbondlist[m].M2].vertype[new_hbondlist[m].arm2]='I';
                    S.M[new_hbondlist[m].M2].nbonds+=1;
                    cout<<"form a hbond"<<endl;
                }
            }
        }
        
        
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
    /*for(int i;i<S.NMOL;i++)
    {
        for(int j=0;j<S.M[i].nbonds;j++)
        {
            totalenergy+=E.hbonde(S.M[i],S.M[S.M[i].hbond_list[j].M2],j);
            //totalenergy+=E.LJ(S.M[i],S.M[j]);
        }
        
    }*/
    return totalenergy;
}
