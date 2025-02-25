//NANOROD: mc.cpp MC Class Function Definitions (Revision Date: October 27,2023)
#include "mc.h" 

void ValidateGrid(const Grid& g, int gridIndex) {
    if(g.n != g.plist.size()) {
        cout << "Grid " << gridIndex << " inconsistency: n=" << g.n 
             << " plist.size()=" << g.plist.size() << endl;
        exit(1);
    }
}

void MC::Sweep()
{
    // Pnew=S.P;
    // E.L=S.L;
    S.WriteDump(0);
    S.WriteData(0);
    double accept=0.0;
   
    int nsample = S.nsweep/1000;
    if(nsample<1)
      nsample=1;
    
    // energy=TotalEnergy();
   
    
    WriteTemplate();
	LogProfile(0,accept);
    
    for(int i=1; i<=S.nsweep; i++)
    {
        time += S.deltat;
        accept += MoveParticle();
        
        if(i%nsample == 0)
        {
            
            LogProfile(i,accept);
            // WriteEnergy(i);
            // S.WriteMol2(i);
            S.WriteDump(i);
            // S.WriteBond(i);
            // S.WriteGrid(i);
            accept=0.0;
        }
    }
	
   
   // S.WriteMol2(S.nsweep);
   S.WriteDump(S.nsweep);
}

void MC::WriteTemplate()
{
    string FileName = "logs/" + S.Description + "_MC.log";
    ofstream out;
    out.open(FileName,ios::app);
    out<<setw(12)<<"sweep"<<"\t"<<setw(12)<<"time"<<"\t"<<setw(12)<<"N_clusters"<<setw(12)<<"Accept"<<"\t"<<endl;
    out.close();
}

void MC::LogProfile(int i, double accept)
{
    string FileName = "logs/" + S.Description + "_MC.log";
    ofstream out;
    out.open(FileName, ios::app);
    out<<setw(12)<<i<<"\t"<<setw(12)<<time<<"\t"<<setw(12)<<S.Ag.size()<<setw(12)<<accept<<"\t"<<endl;
    out.close();
}



//Returns acceptance fraction
double MC::MoveParticle()
{
    

    

    //check if the total particle in the celllists add to N
    
    //check if the aggregate particle number is right
    // int total_particles_ag=0;
    // for(int i=0; i<S.Ag.size(); i++)
    // {
    //     if(S.Ag[i].n != S.Ag[i].plist.size())
    //     {
    //         cout<<"Aggregate particle number wrong"<<endl;
    //         exit(0);
    //     }
    //     total_particles_ag+=S.Ag[i].n;
    // }
    // if(total_particles_ag!=S.NMOL)
    // {
    //     cout<<"Total aggregate particle number wrong"<<endl;
    //     exit(0);
    // }
    // ofstream out;
    // out.open("events.txt",ios::app);
    // ofstream out2;
    // out2.open("Error.txt",ios::app);
 
    double accept=0.0;//accept events
    
    // int N_break;
    // int index;
    // double delta;
    // //break bonds
    // double rand=gsl_rng_uniform(S.gsl_r);
    // if(rand<S.omega_B);//the relative frequency of breaking bond to frequency of diffusion step
    // {
    // for(it=S.H.begin();it!=S.H.end();)
    // {
    //     hbond old_hbond=*it;
    //     //calculate bond_dissociation energy
    //     double E_dis=0;
    //     int free_bonds=0;
    //     E_dis+=S.E_1;//the basic enthalpy change of one bond
    //     //count # of freed bonds
    //     //find neighbor arms,first the one of M1, then the one of M2
    //     int neighborarm1=neighborarm(old_hbond.arm1);
    //     free_bonds+=2;
    //     int bonded_index1;
    //     vector<hbond> old_hbondlist=S.M[old_hbond.M1].hbond_list;
    //     for(int p=0;p<old_hbondlist.size();p++)
    //     {
    //         if(old_hbondlist[p].arm1==neighborarm1)
    //         {
    //             free_bonds-=2;
    //         }
    //         if(old_hbondlist[p].arm1==old_hbond.arm1)
    //         {
    //             bonded_index1=p;
    //         }
    //     } 
    //     int neighborarm2=neighborarm(old_hbond.arm2);
    //     free_bonds+=2;
    //     vector<hbond> bonded_neighbor_hbondlist=S.M[old_hbond.M2].hbond_list;
    //     int bonded_index2;
    //     for(int p=0;p<bonded_neighbor_hbondlist.size();p++)
    //     {
    //         if(bonded_neighbor_hbondlist[p].arm1==neighborarm2)
    //         {
    //             free_bonds-=2;
    //         }
    //         if(bonded_neighbor_hbondlist[p].arm1==old_hbond.arm2)
    //         {
    //             bonded_index2=p;
    //         }
    //     }                 
    //     E_dis+=free_bonds*S.free_bond_freeenergy;
            
    //     if(Arrhenius(1,E_dis,gsl_rng_uniform(S.gsl_r)))
    //     {
    //         //break bond
    //         S.M[old_hbond.M1].vertype[old_hbond.arm1]='A';
    //         out<<"Break bond"<<setw(12)<<old_hbond.M1<<setw(12)<<old_hbond.M2<<setw(12)<<old_hbond.arm1<<setw(12)<<old_hbond.arm2<<endl;
    //         S.M[old_hbond.M1].hbond_list[bonded_index1]=S.M[old_hbond.M1].hbond_list.back();
    //         S.M[old_hbond.M1].nbonds-=1;
            
    //         S.M[old_hbond.M1].hbond_list.pop_back();
    //         S.M[old_hbond.M2].vertype[old_hbond.arm2]='A';
            
    //         S.M[old_hbond.M2].hbond_list[bonded_index2]=S.M[old_hbond.M2].hbond_list.back();
    //         S.M[old_hbond.M2].nbonds-=1;
            
    //         S.M[old_hbond.M2].hbond_list.pop_back();
    //         it=S.H.erase(it);
                
    //     }
    //     else{
    //         ++it;
    //     }
    // }
    // }
    //diffusion rotation
    //diffusion translation
    //diffusion limited bond formation
    int index;
    int N_ag_start = S.Ag.size();
    for(int i = 0; i < N_ag_start; i++)
    {
        int N_ag = S.Ag.size();
        if(N_ag == 1){
            cout<<"Percolated"<<endl;
            break;
        }
        index = gsl_rng_uniform_int(S.gsl_r,N_ag);
        Aggregate new_ag = S.Ag[index];
        double stepsize = S.MCstep/pow(new_ag.n,1.0/3.0);
        XYZ translate_step = RandomTranslate(XYZ(0,0,0),stepsize,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        new_ag.cm = translate_step + new_ag.cm;
        // Create particle list for the new aggregate
        vector<Particle> new_particles;
        for(int j=0; j<new_ag.n; j++) 
        {
            int pid = new_ag.plist[j];
            Particle p = S.P[pid];
            // Update position and grid
            p.position = RandomTranslate(p.position,stepsize,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r)) + translate_step;
            XYZ image_position = image(p.position,S.L);
            p.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.L);
            new_particles.push_back(p);
        }
        // Calculate the energy difference and glauber acceptance
        double delta_energy = 0;
        for(int j=0; j<new_ag.n; j++) {
            int pid = new_ag.plist[j];
            Particle old_particle = S.P[pid];
            int newgID = old_particle.gID;
            int pindex = old_particle.P_ID;
            for(int k=0; k<S.G[newgID].nbr.size(); k++)
            {
                int ngID = S.G[newgID].nbr[k];
                if(S.G[ngID].plist.empty())
                {
                    continue;
                }
                //iterate about molecules in the neighborlist
                list<int>::iterator it;
                for(it=S.G[ngID].plist.begin();it!=S.G[ngID].plist.end();it++)
                {
                    int l=*it;
                    int A_ID1 = S.P[new_particles[j].P_ID].A_ID;
                    int A_ID2 = S.P[l].A_ID;
                    if(A_ID1 == A_ID2)
                    {
                        continue;
                    }
                    double new_r2 = min_d2(new_particles[j].position,S.P[l].position,S.L);
                    double old_r2 = min_d2(old_particle.position,S.P[l].position,S.L);
                    double de = E.Debye_Huckel(new_r2) - E.Debye_Huckel(old_r2);
                    delta_energy += de;

                }
            }
        }
        double rand = gsl_rng_uniform(S.gsl_r);
        if(Glauber(delta_energy,rand)) 
        {
            accept += 1.0;
            S.Ag[index] = new_ag;
            // XYZ image_center = image(new_ag.cm,S.L);
            // Move particles in the aggregate
            
            // First remove particles from their old grids
            for(int j=0; j<new_ag.n; j++) {
                int pid = new_ag.plist[j];
                Particle p = S.P[pid];
                int old_gID = p.gID;
                
                // Remove from old grid first
                S.G[old_gID].n -= 1;
                S.G[old_gID].plist.remove(pid);
            }

            // Then add particles to their new grids
            for(int j=0; j<new_ag.n; j++) 
            {
                S.P[new_particles[j].P_ID] = new_particles[j];
                // Add to new grid
                S.G[new_particles[j].gID].n += 1;
                S.G[new_particles[j].gID].plist.push_back(new_particles[j].P_ID);
            }

            // Search for aggregate in the neighbor of the new particles
            for(int j=0; j<new_ag.n; j++)
            {
                int newgID = new_particles[j].gID;
                int pindex = new_particles[j].P_ID;
                
                
                for(int k=0; k<S.G[newgID].nbr.size(); k++)
                {
                    int ngID = S.G[newgID].nbr[k];
                    // if(ngID < 0 || ngID >= S.NGRID3) 
                    // {
                    //     cout << "Invalid neighbor grid ID: " << ngID << endl;
                    //     exit(1);
                    // }
                    if(S.G[ngID].plist.empty())
                    {
                        continue;
                    }
                    //iterate about molecules in the neighborlist
                    list<int>::iterator it;
                    for(it=S.G[ngID].plist.begin();it!=S.G[ngID].plist.end();it++)
                    {
                        int l=*it;
                        int A_ID1 = S.P[new_particles[j].P_ID].A_ID;
                        int A_ID2 = S.P[l].A_ID;
                        if(A_ID1 == A_ID2)
                        {
                            continue;
                        }
                        double r2 = min_d2(new_particles[j].position,S.P[l].position,S.L);
                        if(r2<S.search2_cm)
                        {
                            // if(A_ID1 < 0 || A_ID1 >= S.Ag.size() ||
                            //    A_ID2 < 0 || A_ID2 >= S.Ag.size()) {
                            //     cout << "Invalid aggregate IDs during merge: " << A_ID1 
                            //          << " " << A_ID2 << " " << S.Ag.size() << endl;
                            //     exit(1);
                            // }

                            // Store the particle list from the old aggregate before any modifications
                            vector<int> old_ag_particles = S.Ag[A_ID2].plist;
                            
                            // Create combined aggregate
                            Aggregate combine_ag = S.Ag[A_ID1];
                            combine_ag.n = S.Ag[A_ID2].n + combine_ag.n;
                            combine_ag.cm = XYZ((S.Ag[A_ID2].cm.x * S.Ag[A_ID2].n + combine_ag.cm.x * combine_ag.n) / combine_ag.n,
                                            (S.Ag[A_ID2].cm.y * S.Ag[A_ID2].n + combine_ag.cm.y * combine_ag.n) / combine_ag.n,
                                            (S.Ag[A_ID2].cm.z * S.Ag[A_ID2].n + combine_ag.cm.z * combine_ag.n) / combine_ag.n);

                            // Update particle list in combined aggregate
                            combine_ag.plist.insert(combine_ag.plist.end(), old_ag_particles.begin(), old_ag_particles.end());
                            
                            // Update A_IDs for all particles from old aggregate
                            for(int pid : old_ag_particles) 
                            {
                                S.P[pid].A_ID = A_ID1;
                            }

                            // Save combined aggregate
                            S.Ag[A_ID1] = combine_ag;

                            // Remove old aggregate
                            S.Ag.erase(S.Ag.begin() + A_ID2);
                            
                            // Update A_IDs for all particles that were in aggregates after A_ID2
                            for(int n=0; n<S.NMOL; n++) 
                            {
                                // if (S.P[n].A_ID == A_ID2)
                                // {
                                //     cout<<"Error: particle "<<n<<" is in the old aggregate "<<A_ID2<<endl;
                                //     exit(1);
                                // }
                                if(S.P[n].A_ID > A_ID2) 
                                {
                                    S.P[n].A_ID -= 1;
                                    // if(S.P[n].A_ID >= S.Ag.size())
                                    // {
                                    //     cout<<"Error: particle "<<n<<" AID exceeds limit "<<S.P[n].A_ID<<endl;
                                    //     exit(1);
                                    // }
                                }   
                            }    
                        }
                    }    
                }  
            }  
        }

        // newmolecule.orientation=RandomRotate(S.M[index].orientation,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        // newmolecule.UpdateVertices();
  
        // int new_hbond=-1;
        // vector<hbond> new_hbondlist;
        // double r2_newbond=S.L*S.L;
        
        // if (newmolecule.nbonds>0)   //calculate hbond_energy
        // {
            
        //     for(int n=0; n<newmolecule.nbonds; n++)
        //     {
        //         double delta_fene=E.hbonde_fene(newmolecule,S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2)-E.hbonde_fene(S.M[index],S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2);
        //         double delta_angle=E.hbonde_angle(newmolecule,S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2)-E.hbonde_angle(S.M[index],S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2);
        //         double delta_dihedral=E.hbonde_dihedral(newmolecule,S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2)-E.hbonde_dihedral(S.M[index],S.M[newmolecule.hbond_list[n].M2],newmolecule.hbond_list[n].arm1,newmolecule.hbond_list[n].arm2);
        //         delta+=delta_fene+delta_angle+delta_dihedral;
        //         //cout<<delta<<endl;
        //     }
        // }
        //add wca energy
        
        //add new wca
    //     for(int l=0;l<nbr_g;l++)
    //     {
    //         //ID of new neighboring grids
    //         int temp_gID=S.G[new_gID].nbr[l];
    //         if(temp_gID>=S.NGRID3)
    //         {
    //             cout<<new_gID<<"Error: wrong grid id"<<endl;
    //             cout<<temp_gID<<endl;
    //             exit(0);
    //         }
    //         if(S.G[temp_gID].plist.empty())
    //         {
    //             continue;
    //         }
    //         //iterate about molecules in the neighborlist
    //         list<int>::iterator it;
    //         for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
    //         {
    //             int j=*it;
    //             if(j!=index)
    //             {
    //                 //image distance
    //                 double rnew2=min_d2(S.M[j].centre,newmolecule.centre,S.L);
    //                 delta+=E.WCA(rnew2);
    //                 //find possible new bonds
    //                 double r2_cm=min_d2(newmolecule.centre,S.M[j].centre,S.L);
    //                 if (r2_cm<S.search2_cm)//check if the cm distance is in the range, can save time
    //                 {
                        
    //                     //iterate about new molecule vertices
    //                     for(int k=0;k<newmolecule.N_VER;k++)
    //                     {
    //                         //iterate about the neighbor vertices
    //                         for(int l=0;l<S.M[j].N_VER;l++)
    //                         {
    //                             double r2_arms=min_d2(newmolecule.ver[k],S.M[j].ver[l],S.L);
    //                             XYZ arm1=image(newmolecule.ver[k],S.L);
    //                             XYZ centre1=image(newmolecule.centre,S.L);
    //                             XYZ neighborarm1=image(newmolecule.ver[neighborarm(k)],S.L);
    //                             XYZ arm2=image(S.M[j].ver[l],S.L);
    //                             XYZ centre2=image(S.M[j].centre,S.L);
    //                             XYZ neighborarm2=image(S.M[j].ver[neighborarm(l)],S.L);
                                
    //                             double alpha=angle_vectors(real_vector(arm1-centre1,S.L),real_vector(arm2-arm1,S.L));
    //                             double beta=angle_vectors(real_vector(arm2-centre2,S.L),real_vector(arm1-arm2,S.L));
    //                             double xhi=dihedral_vectors(real_vector(neighborarm1-centre1,S.L),real_vector(centre2-centre1,S.L),real_vector(neighborarm2-centre2,S.L));
                              
    //                             //check type, distance between arms, angles and dihedral
    //                             if (newmolecule.vertype[k]=='A'&&S.M[j].vertype[l]=='A'&&r2_arms<S.sfearchl2_bond&&alpha<S.search_angle&&beta<S.search_angle&& abs(xhi)<S.search_xhi)
    //                             {
                                    
                                    
    //                                 bool bonded=false;//check if there is already bonds between the two molecules, only one bond can exist between two molecules
    //                                 for(int n=0;n<old_hbondlist.size();n++)
    //                                 {
    //                                     if(old_hbondlist[n].M2==j)
    //                                     {
    //                                         bonded=true;
    //                                     }
    //                                 }
    //                                 for(int n=0;n<new_hbondlist.size();n++)
    //                                 {
    //                                     if(new_hbondlist[n].M2==j)
    //                                     {
    //                                         bonded=true;
    //                                     }
    //                                 }
    //                                 if(bonded==false)
    //                                 {
    //                                     new_hbondlist.push_back(hbond(index,j,k,l));
    //                                     new_hbond=1;
    //                                 }
    //                             }
    //                         }
    //                     }
                        
    //                 }


    //             }
    //         }
    //     }
    //     //subtract old wca
    //     int old_gID=S.M[index].gID;
    //     for(int l=0;l<nbr_g;l++)
    //     {
    //         //ID of new neighboring grids
    //         int temp_gID=S.G[old_gID].nbr[l];
    //         list<int>::iterator it;
    //         for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
    //         {
    //             int j=*it;
    //             if(j!=index)
    //             {
    //                 //image distance
    //                 double rold2=min_d2(S.M[j].centre,S.M[index].centre,S.L);
    //                 delta-=E.WCA(rold2);
                    


    //             }
    //         }
    //     }
    //     // if(Glauber(delta,gsl_rng_uniform(S.gsl_r)))
    //     // {
    //         //Accept move
    //         S.M[index]=newmolecule;
    
    //         accept+=1.0;
    //         energy+=delta;

    //         //update grid list
    //         if(new_gID!=old_gID)
    //         {
    //             S.G[old_gID].n-=1;
    //             S.G[old_gID].plist.remove(S.M[index].MOL_ID);
                
    //             S.G[new_gID].n+=1;
    //             S.G[new_gID].plist.push_back(S.M[index].MOL_ID);
    //         }
        
    //         //form bonds if accept
    //         if (new_hbond!=-1)
    //         {
    //             for(int m=0;m<new_hbondlist.size();m++) 
    //             {
    //                 S.M[index].hbond_list.push_back(new_hbondlist[m]);
    //                 S.M[index].vertype[new_hbondlist[m].arm1]='I';
    //                 S.M[index].nbonds+=1;
    //                 S.M[new_hbondlist[m].M2].hbond_list.push_back(hbond(new_hbondlist[m].M2,new_hbondlist[m].M1,new_hbondlist[m].arm2,new_hbondlist[m].arm1));
    //                 S.M[new_hbondlist[m].M2].vertype[new_hbondlist[m].arm2]='I';
    //                 S.M[new_hbondlist[m].M2].nbonds+=1;
    //                 out<<"form a hbond"<<setw(12)<<new_hbondlist[m].M1<<setw(12)<<new_hbondlist[m].M2<<setw(12)<<new_hbondlist[m].arm1<<setw(12)<<new_hbondlist[m].arm2<<endl;
    //                 if(new_hbondlist[m].M1<new_hbondlist[m].M2)
    //                 {
    //                     S.H.push_back(new_hbondlist[m]);
    //                 }
    //                 else
    //                 {
    //                     S.H.push_back(hbond(new_hbondlist[m].M2,new_hbondlist[m].M1,new_hbondlist[m].arm2,new_hbondlist[m].arm1));
    //                 }
    //             }

    //         }
    //     // }
        
        
    }
    for(int i=0; i<S.NMOL; i++)
    {
        index = gsl_rng_uniform_int(S.gsl_r,S.NMOL);
        Particle new_particle = S.P[index];
        new_particle.position = RandomTranslate(new_particle.position,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        XYZ image_position = image(new_particle.position,S.L);
        new_particle.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.L);
        // Calculate the energy difference and glauber acceptance
        double delta_energy = 0;
        Particle old_particle = S.P[index];
        int oldgID = old_particle.gID;
        for(int k=0; k<S.G[oldgID].nbr.size(); k++)
        {
            int ngID = S.G[oldgID].nbr[k];
            if(S.G[ngID].plist.empty())
            {
                continue;
            }
            //iterate about molecules in the neighborlist
            list<int>::iterator it;
            for(it=S.G[ngID].plist.begin();it!=S.G[ngID].plist.end();it++)
            {
                int l=*it;
                double new_r2 = min_d2(new_particle.position,S.P[l].position,S.L);
                double old_r2 = min_d2(old_particle.position,S.P[l].position,S.L);
                double de = E.total_energy(new_r2) - E.total_energy(old_r2);
                delta_energy += de;

            }
        }
        if(Glauber(delta_energy,gsl_rng_uniform(S.gsl_r)))
        {
            accept += 1.0;
            energy += delta_energy;
            // First remove particles from their old grids
            
           
                
                
            // Remove from old grid first
            S.G[oldgID].n -= 1;
            S.G[oldgID].plist.remove(index);
        

            // Then add particles to their new grids
            
            S.P[index] = new_particle;
            // Add to new grid
            S.G[new_particle.gID].n += 1;
            S.G[new_particle.gID].plist.push_back(index);
  
        }
    }
    // for(int i=0; i<S.NGRID3; i++) {
    //     ValidateGrid(S.G[i], i);
    // }

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
// bool MC::Arrhenius(double A,double delta, double rand)
// {
    
//     if(A*(exp(-delta))>rand)
//         {
//         return true;}
//     else
//         return false;
    
    
// }
// double MC::WCAEnergy()
// {
//     double ewca=0;
//     list<int>::iterator it;
//     for(int i=0;i<S.NMOL;i++)
//     {
//         int old_gID=S.M[i].gID;
//         for(int l=0;l<nbr_g;l++)
//         {
//             //ID of new neighboring grids
//             int temp_gID=S.G[old_gID].nbr[l];
            
//             for(it=S.G[temp_gID].plist.begin();it!=S.G[temp_gID].plist.end();it++)
//             {
//                 int j=*it;
//                 if(j>i)//avoid double counting
//                 {
//                     //image distance
//                     double r2=min_d2(S.M[j].centre,S.M[i].centre,S.L);
//                     ewca+=E.WCA(r2);
                    


//                 }
//             }
//         }
//     }
//     return ewca;
// }
// // double MC::FENE_energy()
// {
//     double efene=0;
//     list<hbond>::iterator it;
//     for(it=S.H.begin();it!=S.H.end();it++)
//     {
//         hbond j=*it;
//         efene+=E.hbonde_fene(S.M[j.M1],S.M[j.M2],j.arm1,j.arm2);
//     }
//     return efene;
// }
// double MC::Angle_energy()
// {
//     double eangle=0;
//     list<hbond>::iterator it;
//     for(it=S.H.begin();it!=S.H.end();it++)
//     {
//         hbond j=*it;
//         eangle+=E.hbonde_angle(S.M[j.M1],S.M[j.M2],j.arm1,j.arm2);
//     }
//     return eangle;
// }
// double MC::Dihedral_energy()
// {
//     double edihedral=0;
//     list<hbond>::iterator it;
//     for(it=S.H.begin();it!=S.H.end();it++)
//     {
//         hbond j=*it;
//         edihedral+=E.hbonde_dihedral(S.M[j.M1],S.M[j.M2],j.arm1,j.arm2);
//     }
//     return edihedral;   
// }
// double MC::bond_energy()
// {
//     double ebond;
//     ebond=-S.H.size()*S.E_1;
//     return ebond;
// }
// double MC::bond_freeze_freenerngy()
// {
//     double efreeze;
//     int nfreeze;
//     for(int i=0;i<S.NMOL;i++)
//     {
//         Molecule m=S.M[i];
//         for(int l=0;l<m.N_VER;l++)
//         {
//             if(m.vertype[l]=='I'||m.vertype[neighborarm(l)]=='I')
//             nfreeze+=1;

//         }
//     }
//     efreeze=nfreeze*S.free_bond_freeenergy;
// }
// double MC::TotalEnergy()
// {
//     double totalenergy=0.0;
//     /*for(int i;i<S.NMOL;i++)
//     {
//         for(int j=0;j<S.M[i].nbonds;j++)
//         {
//             totalenergy+=E.hbonde(S.M[i],S.M[S.M[i].hbond_list[j].M2],j);
//             //totalenergy+=E.LJ(S.M[i],S.M[j]);
//         }
        
//     }*/
//     return totalenergy;
// }
// double MC::WriteEnergy(int timestep)
// {
//     ofstream out;
//     out.open("energy.txt",ios::app);
//     out<<"TIMESTEP"<<endl;
//     out<<timestep<<endl;
//     out<<"WCA"<<'\t'<<"FENE"<<'\t'<<"Angle"<<'\t'<<"Dihedral"<<'\t'<<"Bond"<<'\t'<<"Bond_freeze"<<endl;
//     out<<WCAEnergy()<<'\t'<<FENE_energy()<<'\t'<<Angle_energy()<<'\t'<<Dihedral_energy()<<'\t'<<bond_energy()<<'\t'<<bond_freeze_freenerngy()<<endl;
//     out.close();
// }
