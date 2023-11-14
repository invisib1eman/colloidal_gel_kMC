//NANOROD: mc.cpp MC Class Function Definitions (Revision Date: October 27,2023)
#include "mc.h" 
void MC::Sweep()
{
    Mnew=S.M;
    E.L=S.L;
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
    out<<"bondlist"<<endl;
    for(int j=0;j<S.NMOL;j++)
    {
        if(S.M[j].hbond_list.size()>0)
        {
            for(int k=0;k<S.M[j].hbond_list.size();k++)
            {
                out<<setw(12)<<S.M[j].hbond_list[k].M1<<setw(12)<<S.M[j].hbond_list[k].M2<<setw(12)<<S.M[j].hbond_list[k].arm1<<setw(12)<<S.M[j].hbond_list[k].arm2<<endl;
            }
        }
        
    }
    out.close();
}


//Returns acceptance fraction
double MC::MoveMolecule()
{
    
    double accept=0.0;//accept events
    list<int>::iterator it;
    
    //vector<int> sequence_M=generateRandom(S.NMOL);//randomly generate a sequence of the N molecules;
    for(int i=0; i<S.NMOL; i++)
    {
        int index=gsl_rng_uniform_int(S.gsl_r,S.NMOL);
        
        double delta=0;//energy difference
        //int index=sequence_M[i];//index of the current trial molecule
        Molecule newmolecule=S.M[index];
        //Check the molecule bonds and see if the bonds will break, the probability of breaking a bond is given by Arrhenius formula
        vector<hbond> old_hbondlist=S.M[index].hbond_list;
        for(int n=0;n<old_hbondlist.size();n++)
        {
            hbond old_hbond=old_hbondlist[n];
            //calculate bond_dissociation energy
            double E_dis;
            int free_bonds=0;
            E_dis+=S.E_1;//the basic enthalpy change of one bond
            //count # of freed bonds
            //find neighbor arms,first the one of M1, then the one of M2
            int neighborarm1=neighborarm(old_hbond.arm1);
            free_bonds+=2;
            for(int p=0;p<old_hbondlist.size();p++)
            {
                if(old_hbondlist[p].arm1==neighborarm1)
                {
                    free_bonds-=2;
                }
            } 
            int neighborarm2=neighborarm(old_hbond.arm2);
            free_bonds+=2;
            vector<hbond> bonded_neighbor_hbondlist=S.M[old_hbond.M2].hbond_list;
            for(int p=0;p<bonded_neighbor_hbondlist.size();p++)
            {
                if(bonded_neighbor_hbondlist[p].arm1==neighborarm2)
                {
                    free_bonds-=2;
                }
            }                 
            E_dis+=free_bonds*S.free_bond_freeenergy;

            
        }
        newmolecule.centre=RandomTranslate(S.M[index].centre,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        XYZ image_center=image(newmolecule.centre,S.L);
        int new_gID=GridIndex_xyz(image_center,S.NGRID,S.GRIDL,S.L);
        newmolecule.gID=new_gID;

        newmolecule.orientation=RandomRotate(S.M[index].orientation,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        newmolecule.UpdateVertices();
        int new_hbond=-1;
        vector<hbond> new_hbondlist;
        double r2_newbond=S.L*S.L;
        if (newmolecule.nbonds>0)   //calculate hbond_energy
        {
            
            for(int n=0; n<newmolecule.nbonds; n++)
            {
                delta+=E.hbonde_fene(newmolecule,S.M[newmolecule.hbond_list[n].M2],n)-E.hbonde_fene(S.M[index],S.M[newmolecule.hbond_list[n].M2],n);
                //cout<<delta<<endl;
            }
        }
        //add wca energy
        
        //add new wca
        for(int l=0;l<nbr_g;l++)
        {
            //ID of new neighboring grids
            int temp_gID=S.G[new_gID].nbr[l];
            for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
            {
                int j=*it;
                if(j!=index)
                {
                    //image distance
                    double rnew2=min_d2(S.M[j].centre,newmolecule.centre,S.L);
                    delta+=E.WCA(rnew2);
                    //find possible new bonds
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
                                    
                                    bool bonded=false;//check if there is already bonds between the two molecules, only one bond can exist between two molecules
                                    for(int n=0;n<old_hbondlist.size();n++)
                                    {
                                        if(old_hbondlist[n].M2==j)
                                        {
                                            bonded=true;
                                        }
                                    }
                                    for(int n=0;n<new_hbondlist.size();n++)
                                    {
                                        if(new_hbondlist[n].M2==j)
                                        {
                                            bonded=true;
                                        }
                                    }
                                    if(bonded==false)
                                    {
                                        new_hbondlist.push_back(hbond(index,j,k,l));
                                        new_hbond=1;
                                    }
                                }
                            }
                        }
                        
                    }


                }
            }
        }
        //subtract old wca
        int old_gID=S.M[index].gID;
        for(int l=0;l<nbr_g;l++)
        {
            //ID of new neighboring grids
            int temp_gID=S.G[old_gID].nbr[l];
            for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
            {
                int j=*it;
                if(j!=index)
                {
                    //image distance
                    double rold2=min_d2(S.M[j].centre,S.M[index].centre,S.L);
                    delta-=E.WCA(rold2);
                    


                }
            }
        }
        if(Glauber(delta,gsl_rng_uniform(S.gsl_r)))
        {
            //Accept move
            S.M[index]=newmolecule;
            cout<<"accept"<<endl;
            accept+=1.0;
            energy+=delta;

            //update grid list
            if(new_gID!=old_gID)
            {
                S.G[old_gID].n-=1;
                S.G[old_gID].plist.remove(S.M[index].MOL_ID);
                cout<<index<<S.M[index].MOL_ID;
                S.G[new_gID].n+=1;
                S.G[new_gID].plist.push_back(S.M[index].MOL_ID);
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
bool MC::Arrhenius(double delta, double rand)
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
