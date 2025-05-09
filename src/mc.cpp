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
    S.WriteDump(0);
    S.WriteData(0);
    double accept=0.0;
    Find_Neighbors();
    Cluster_Particles();
    // Output 1000 frames
    int nsample = S.nsweep/N_frame;
    if(nsample<1)
      nsample=1;
    
    WriteTemplate();
	LogProfile(0,accept);
    // Sweep
    if (S.mode == "single_particle")
    {
        for(int i=1; i<=S.nsweep; i++)
        {
            time += S.deltat; 
            accept += MoveParticle_Single_Particle();
            if(i%nsample == 0)
            {
                
                LogProfile(i,accept);
                S.WriteDump(i);
                accept=0.0;
            }
        }
    }
    if (S.mode == "cluster_free_roll")
    {
        for(int i=1; i<=S.nsweep; i++)
        {
            time += S.deltat; 
            accept += MoveParticle_Cluster_Free_Roll();
            
            if(i%nsample == 0)
            {
                
                LogProfile(i,accept);
                S.WriteDump(i);
                accept=0.0;
            }
        }   
    }
    if (S.mode == "cluster_rigid")
    {
        for(int i=1; i<=S.nsweep; i++)
        {
            time += S.deltat; 
            accept += MoveParticle_Cluster_Rigid();
            
            if(i%nsample == 0)
            {
                
                LogProfile(i,accept);
                S.WriteDump(i);
                accept=0.0;
            }
        }
    }
    if (S.mode == "cluster_alpha")
    {
        for(int i=1; i<=S.nsweep; i++)
        {
            time += S.deltat; 
            accept += MoveParticle_Cluster_Alpha();
            if(i%nsample == 0)
            {
                
                LogProfile(i,accept);
                S.WriteDump(i);
                accept=0.0;
            }
        }
    }
    if (S.mode == "cluster_alpha_morse")
    {
        for(int i=1; i<=S.nsweep; i++)
        {
            time += S.deltat; 
            accept += MoveParticle_Cluster_Alpha_Morse();
            if(i%nsample == 0)
            {
                
                LogProfile(i,accept);
                S.WriteDump(i);
                accept=0.0;
            }
        }
    }
    S.WriteDump(S.nsweep);
}

void MC::WriteTemplate()
{
    string FileName = S.log_file_name;
    ofstream out;
    out.open(FileName,ios::trunc);
    out<<setw(12)<<"sweep"<<"\t"<<setw(12)<<"time"<<"\t"<<setw(12)<<"N_clusters"<<setw(12)<<"Accept"<<"\t"<<endl;
    out.close();
}

void MC::LogProfile(int i, double accept)
{
    string FileName = S.log_file_name;
    ofstream out;
    out.open(FileName, ios::app);
    out<<setw(12)<<i<<"\t"<<setw(12)<<time<<"\t"<<setw(12)<<S.Ag.size()<<setw(12)<<accept<<"\t"<<endl;
    out.close();
}
//mode: single_particle
double MC::MoveParticle_Single_Particle()
{
    double accept=0.0;//accept events
    int index;
    // Loop over all particles
    for(int i=0; i<S.NMOL; i++)
    {
        index = gsl_rng_uniform_int(S.gsl_r,S.NMOL);
        Particle new_particle = S.P[index];
        new_particle.position = RandomTranslate(new_particle.position,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        XYZ image_position = image(new_particle.position,S.BoxLength);
        new_particle.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.BoxLength);
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
                if(l == index)
                {
                    continue;
                }
                double new_r2 = min_d2(new_particle.position,S.P[l].position,S.BoxLength);
                double old_r2 = min_d2(old_particle.position,S.P[l].position,S.BoxLength);
                if (new_r2 < 1)
                {
                    cout << "Too close" << endl;
                    exit(1);
                }
                double de = E.total_energy(new_r2) - E.total_energy(old_r2);
                delta_energy += de;

            }
        }
        if(Glauber(delta_energy,gsl_rng_uniform(S.gsl_r)))
        {
            accept += 1.0;
            energy += delta_energy;
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
    
    return accept/double(S.NMOL);
}

//mode: cluster_rigid
double MC::MoveParticle_Cluster_Rigid()
{
    double accept=0.0;//accept events
    int index;
    int N_ag_start = S.Ag.size();
    // Loop over all aggregates
    for(int i = 0; i < N_ag_start; i++)
    {
        int N_ag = S.Ag.size();
        // if(N_ag == 1){
        //     cout<<"Percolated"<<endl;
        //     break;
        // }
        index = gsl_rng_uniform_int(S.gsl_r,N_ag);
        Aggregate new_ag = S.Ag[index];
        double stepsize;
        if (S.fake_acceleration == 0)
        {
            stepsize = S.MCstep/pow(new_ag.n,1.0/6.0);
        }
        else
        {
            stepsize = S.MCstep;
        }
        XYZ translate_step = RandomTranslate(XYZ(0,0,0),stepsize,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
        new_ag.cm = translate_step + new_ag.cm;
        // Create particle list for the new aggregate
        vector<Particle> new_particles;
        for(int j=0; j<new_ag.n; j++) 
        {
            int pid = new_ag.plist[j];
            Particle p = S.P[pid];
            // Update position and grid
            p.position = p.position + translate_step;
            XYZ image_position = image(p.position,S.BoxLength);
            p.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.BoxLength);
            new_particles.push_back(p);
        }
        // Calculate the energy difference and glauber acceptance
        double delta_energy = 0;
        for(int j=0; j<new_ag.n; j++) {
                int pID = new_ag.plist[j];
                Particle old_particle = S.P[pID];
                Particle new_particle = new_particles[j];
                int oldgID = old_particle.gID;
                int newgID = new_particle.gID;
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
                        int A_ID1 = S.P[pID].A_ID;
                        int A_ID2 = S.P[l].A_ID;
                        if(A_ID1 == A_ID2)
                        {
                            continue;
                        }
                        double old_r2 = min_d2(old_particle.position,S.P[l].position,S.BoxLength);
                        double de = - E.total_energy(old_r2);
                        delta_energy += de;
                    }
                }
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
                        int A_ID1 = S.P[pID].A_ID;
                        int A_ID2 = S.P[l].A_ID;
                        if(A_ID1 == A_ID2)
                        {
                            continue;
                        }
                        double new_r2 = min_d2(new_particles[j].position,S.P[l].position,S.BoxLength);
                        double de = E.total_energy(new_r2);
                        delta_energy += de;
                    }
                }
            }
        double rand = gsl_rng_uniform(S.gsl_r);
        if(Glauber(delta_energy,rand)) 
        {
            accept += 1.0;
            S.Ag[index] = new_ag;
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
                        double r2 = min_d2(new_particles[j].position,S.P[l].position,S.BoxLength);
                        if(r2<S.search2_cm)
                        {
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
                                if(S.P[n].A_ID > A_ID2) 
                                {
                                    S.P[n].A_ID -= 1;
                                }   
                            }    
                        }
                    }    
                }  
            }  
        }   
    }
    return accept/double(S.NMOL);
}
//mode: cluster_free_roll
double MC::MoveParticle_Cluster_Free_Roll()
{
    double accept=0.0;//accept events
    int index;
    int N_ag_start = S.Ag.size();
    double rand = gsl_rng_uniform(S.gsl_r);
    // Determine if the move is cluster diffuse or single particle crawl (relaxation)
    // Aggregate diffusion happens with probability 1 - p_crawl
    if(rand > S.p_crawl)
    {
        // Loop over all aggregates
        for(int i = 0; i < N_ag_start; i++)
        {
            int N_ag = S.Ag.size();
            // if(N_ag == 1){
            //     cout<<"Percolated"<<endl;
            //     break;
            // }
            index = gsl_rng_uniform_int(S.gsl_r,N_ag);
            Aggregate new_ag = S.Ag[index];
            double stepsize;
            if (S.fake_acceleration == 0)
            {
                stepsize = S.MCstep/pow(new_ag.n,1.0/6.0);
            }
            else
            {
                stepsize = S.MCstep;
            }
            XYZ translate_step = RandomTranslate(XYZ(0,0,0),stepsize,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
            new_ag.cm = translate_step + new_ag.cm;
            // Create particle list for the new aggregate
            vector<Particle> new_particles;
            for(int j=0; j<new_ag.n; j++) 
            {
                int pid = new_ag.plist[j];
                Particle p = S.P[pid];
                // Update position and grid
                p.position = p.position + translate_step;
                XYZ image_position = image(p.position,S.BoxLength);
                p.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.BoxLength);
                new_particles.push_back(p);
            }
            // Calculate the energy difference and glauber acceptance
            double delta_energy = 0;
            for(int j=0; j<new_ag.n; j++) {
                int pID = new_ag.plist[j];
                Particle old_particle = S.P[pID];
                Particle new_particle = new_particles[j];
                int oldgID = old_particle.gID;
                int newgID = new_particle.gID;
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
                        int A_ID1 = S.P[pID].A_ID;
                        int A_ID2 = S.P[l].A_ID;
                        if(A_ID1 == A_ID2)
                        {
                            continue;
                        }
                        double old_r2 = min_d2(old_particle.position,S.P[l].position,S.BoxLength);
                        double de = - E.total_energy(old_r2);
                        delta_energy += de;
                    }
                }
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
                        int A_ID1 = S.P[pID].A_ID;
                        int A_ID2 = S.P[l].A_ID;
                        if(A_ID1 == A_ID2)
                        {
                            continue;
                        }
                        double new_r2 = min_d2(new_particles[j].position,S.P[l].position,S.BoxLength);
                        double de = E.total_energy(new_r2);
                        delta_energy += de;
                    }
                }
            }
            double rand = gsl_rng_uniform(S.gsl_r);
            if(Glauber(delta_energy,rand)) 
            {
                accept += 1.0;
                S.Ag[index] = new_ag;
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
                            double r2 = min_d2(new_particles[j].position,S.P[l].position,S.BoxLength);
                            if(r2<S.search2_cm)
                            {
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
                                    if(S.P[n].A_ID > A_ID2) 
                                    {
                                        S.P[n].A_ID -= 1;
                                    }   
                                }    
                            }
                        }    
                    }  
                }  
            }   
        }
    }
    // Single particle crawl (relaxation) happens with probability p_crawl
    else
    {
        // Loop over all particles
        for(int i=0; i<S.NMOL; i++)
        {
            index = gsl_rng_uniform_int(S.gsl_r,S.NMOL);
            Particle new_particle = S.P[index];
            if (S.Ag[new_particle.A_ID].n == 1) 
            {
                continue;
            }
            new_particle.position = RandomTranslate(new_particle.position,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
            XYZ image_position = image(new_particle.position,S.BoxLength);
            new_particle.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.BoxLength);
            // Calculate the energy difference and glauber acceptance
            double delta_energy = 0;
            Particle old_particle = S.P[index];
            int oldgID = old_particle.gID;
            int newgID = new_particle.gID;
            // Calculate energy difference for the particle move
            // First, remove the current particle from neighbors' nbr_particles lists
            
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
                    if(l == index)
                    {
                        continue;
                    }
                    double old_r2 = min_d2(old_particle.position,S.P[l].position,S.BoxLength);
                    // if (new_r2 < 1)
                    // {
                    //     cout << "Too close" << endl;
                    //     exit(1);
                    // }
                    double de = - E.total_energy(old_r2);
                    delta_energy += de;

                }
            }
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
                    if(l == index)
                    {
                        continue;
                    }
                    double new_r2 = min_d2(new_particle.position,S.P[l].position,S.BoxLength);
                    double de = E.total_energy(new_r2);
                    delta_energy += de;

                }
            }
            if(Glauber(delta_energy,gsl_rng_uniform(S.gsl_r)))
            {
                accept += 1.0;
                energy += delta_energy;    
                // Remove from old grid first
                S.G[oldgID].n -= 1;
                S.G[oldgID].plist.remove(index);
                // Add to new grid
                S.G[new_particle.gID].n += 1;
                S.G[new_particle.gID].plist.push_back(index);
                S.P[index] = new_particle;
    
            }
            
        }
        Find_Neighbors();
        Cluster_Particles();
    }
    return accept/double(S.NMOL);
}

double MC::MoveParticle_Cluster_Alpha()
{
    double accept=0.0;//accept events
    int index;
    int N_ag_start = S.Ag.size();
    double rand = gsl_rng_uniform(S.gsl_r);
    // Determine if the move is cluster diffuse or single particle crawl (relaxation)
    // Aggregate diffusion happens with probability 1 - p_crawl
    if(rand > S.p_crawl)
    {
        // Loop over all aggregates
        for(int i = 0; i < N_ag_start; i++)
        {
            int N_ag = S.Ag.size();
            // if(N_ag == 1){
            //     cout<<"Percolated"<<endl;
            //     break;
            // }
            index = gsl_rng_uniform_int(S.gsl_r,N_ag);
            Aggregate new_ag = S.Ag[index];
            double stepsize;
            if (S.fake_acceleration == 0)
            {
                double exponent_stepsize = S.alpha/2.0;
                stepsize = S.MCstep/pow(new_ag.n,exponent_stepsize);
            }
            else
            {
                stepsize = S.MCstep;
            }
            XYZ translate_step = RandomTranslate(XYZ(0,0,0),stepsize,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
            new_ag.cm = translate_step + new_ag.cm;
            // Create particle list for the new aggregate
            vector<Particle> new_particles;
            for(int j=0; j<new_ag.n; j++) 
            {
                int pid = new_ag.plist[j];
                Particle p = S.P[pid];
                // Update position and grid
                p.position = p.position + translate_step;
                XYZ image_position = image(p.position,S.BoxLength);
                p.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.BoxLength);
                new_particles.push_back(p);
            }
            // Calculate the energy difference and glauber acceptance
            double delta_energy = 0;
            for(int j=0; j<new_ag.n; j++) {
                int pID = new_ag.plist[j];
                Particle old_particle = S.P[pID];
                Particle new_particle = new_particles[j];
                int oldgID = old_particle.gID;
                int newgID = new_particle.gID;
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
                        int A_ID1 = S.P[pID].A_ID;
                        int A_ID2 = S.P[l].A_ID;
                        if(A_ID1 == A_ID2)
                        {
                            continue;
                        }
                        double old_r2 = min_d2(old_particle.position,S.P[l].position,S.BoxLength);
                        double de = - E.total_energy(old_r2);
                        delta_energy += de;
                    }
                }
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
                        int A_ID1 = S.P[pID].A_ID;
                        int A_ID2 = S.P[l].A_ID;
                        if(A_ID1 == A_ID2)
                        {
                            continue;
                        }
                        double new_r2 = min_d2(new_particles[j].position,S.P[l].position,S.BoxLength);
                        double de = E.total_energy(new_r2);
                        delta_energy += de;
                    }
                }
            }
            double rand = gsl_rng_uniform(S.gsl_r);
            if(Glauber(delta_energy,rand)) 
            {
                accept += 1.0;
                S.Ag[index] = new_ag;
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
                            double r2 = min_d2(new_particles[j].position,S.P[l].position,S.BoxLength);
                            if(r2<S.search2_cm)
                            {
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
                                    if(S.P[n].A_ID > A_ID2) 
                                    {
                                        S.P[n].A_ID -= 1;
                                    }   
                                }    
                            }
                        }    
                    }  
                }  
            }   
        }
    }
    // Single particle crawl (relaxation) happens with probability p_crawl
    else
    {
        // Loop over all particles
        for(int i=0; i<S.NMOL; i++)
        {
            index = gsl_rng_uniform_int(S.gsl_r,S.NMOL);
            Particle new_particle = S.P[index];
            if (S.Ag[new_particle.A_ID].n == 1) 
            {
                continue;
            }
            new_particle.position = RandomTranslate(new_particle.position,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
            XYZ image_position = image(new_particle.position,S.BoxLength);
            new_particle.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.BoxLength);
            // Calculate the energy difference and glauber acceptance
            double delta_energy = 0;
            Particle old_particle = S.P[index];
            int oldgID = old_particle.gID;
            int newgID = new_particle.gID;
            // Calculate energy difference for the particle move
            // First, remove the current particle from neighbors' nbr_particles lists
            
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
                    if(l == index)
                    {
                        continue;
                    }
                    double old_r2 = min_d2(old_particle.position,S.P[l].position,S.BoxLength);
                    // if (new_r2 < 1)
                    // {
                    //     cout << "Too close" << endl;
                    //     exit(1);
                    // }
                    double de = - E.total_energy(old_r2);
                    delta_energy += de;

                }
            }
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
                    if(l == index)
                    {
                        continue;
                    }
                    double new_r2 = min_d2(new_particle.position,S.P[l].position,S.BoxLength);
                    double de = E.total_energy(new_r2);
                    delta_energy += de;

                }
            }
            if(Glauber(delta_energy,gsl_rng_uniform(S.gsl_r)))
            {
                accept += 1.0;
                energy += delta_energy;    
                // Remove from old grid first
                S.G[oldgID].n -= 1;
                S.G[oldgID].plist.remove(index);
                // Add to new grid
                S.G[new_particle.gID].n += 1;
                S.G[new_particle.gID].plist.push_back(index);
                S.P[index] = new_particle;
    
            }
            
        }
        Find_Neighbors();
        Cluster_Particles();
    }
    return accept/double(S.NMOL);
}
double MC::MoveParticle_Cluster_Alpha_Morse()
{
    double accept=0.0;//accept events
    int index;
    int N_ag_start = S.Ag.size();
    double rand = gsl_rng_uniform(S.gsl_r);
    // Determine if the move is cluster diffuse or single particle crawl (relaxation)
    // Aggregate diffusion happens with probability 1 - p_crawl
    if(rand > S.p_crawl)
    {
        // Loop over all aggregates
        for(int i = 0; i < N_ag_start; i++)
        {
            int N_ag = S.Ag.size();
            // if(N_ag == 1){
            //     cout<<"Percolated"<<endl;
            //     break;
            // }
            index = gsl_rng_uniform_int(S.gsl_r,N_ag);
            Aggregate new_ag = S.Ag[index];
            double stepsize;
            if (S.fake_acceleration == 0)
            {
                double exponent_stepsize = S.alpha/2.0;
                stepsize = S.MCstep/pow(new_ag.n,exponent_stepsize);
            }
            else
            {
                stepsize = S.MCstep;
            }
            XYZ translate_step = RandomTranslate(XYZ(0,0,0),stepsize,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
            new_ag.cm = translate_step + new_ag.cm;
            // Create particle list for the new aggregate
            vector<Particle> new_particles;
            for(int j=0; j<new_ag.n; j++) 
            {
                int pid = new_ag.plist[j];
                Particle p = S.P[pid];
                // Update position and grid
                p.position = p.position + translate_step;
                XYZ image_position = image(p.position,S.BoxLength);
                p.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.BoxLength);
                new_particles.push_back(p);
            }
            // Calculate the energy difference and glauber acceptance
            double delta_energy = 0;
            for(int j=0; j<new_ag.n; j++) {
                int pID = new_ag.plist[j];
                Particle old_particle = S.P[pID];
                Particle new_particle = new_particles[j];
                int oldgID = old_particle.gID;
                int newgID = new_particle.gID;
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
                        int A_ID1 = S.P[pID].A_ID;
                        int A_ID2 = S.P[l].A_ID;
                        if(A_ID1 == A_ID2)
                        {
                            continue;
                        }
                        double old_r2 = min_d2(old_particle.position,S.P[l].position,S.BoxLength);
                        double de = - E.morse_potential(old_r2)-E.Debye_Huckel(old_r2);
                        delta_energy += de;
                    }
                }
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
                        int A_ID1 = S.P[pID].A_ID;
                        int A_ID2 = S.P[l].A_ID;
                        if(A_ID1 == A_ID2)
                        {
                            continue;
                        }
                        double new_r2 = min_d2(new_particles[j].position,S.P[l].position,S.BoxLength);
                        double de = E.morse_potential(new_r2)+E.Debye_Huckel(new_r2);
                        delta_energy += de;
                    }
                }
            }
            double rand = gsl_rng_uniform(S.gsl_r);
            if(Glauber(delta_energy,rand)) 
            {
                accept += 1.0;
                S.Ag[index] = new_ag;
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
                            double r2 = min_d2(new_particles[j].position,S.P[l].position,S.BoxLength);
                            if(r2<S.search2_cm)
                            {
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
                                    if(S.P[n].A_ID > A_ID2) 
                                    {
                                        S.P[n].A_ID -= 1;
                                    }   
                                }    
                            }
                        }    
                    }  
                }  
            }   
        }
    }
    // Single particle crawl (relaxation) happens with probability p_crawl
    else
    {
        // Loop over all particles
        for(int i=0; i<S.NMOL; i++)
        {
            index = gsl_rng_uniform_int(S.gsl_r,S.NMOL);
            Particle new_particle = S.P[index];
            if (S.Ag[new_particle.A_ID].n == 1) 
            {
                continue;
            }
            new_particle.position = RandomTranslate(new_particle.position,S.MCstep,gsl_rng_uniform(S.gsl_r),gsl_rng_uniform(S.gsl_r));
            XYZ image_position = image(new_particle.position,S.BoxLength);
            new_particle.gID = GridIndex_xyz(image_position,S.NGRID,S.GRIDL,S.BoxLength);
            // Calculate the energy difference and glauber acceptance
            double delta_energy = 0;
            Particle old_particle = S.P[index];
            int oldgID = old_particle.gID;
            int newgID = new_particle.gID;
            // Calculate energy difference for the particle move
            // First, remove the current particle from neighbors' nbr_particles lists
            
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
                    if(l == index)
                    {
                        continue;
                    }
                    double old_r2 = min_d2(old_particle.position,S.P[l].position,S.BoxLength);
                    // if (new_r2 < 1)
                    // {
                    //     cout << "Too close" << endl;
                    //     exit(1);
                    // }
                    double de = - E.morse_potential(old_r2)-E.Debye_Huckel(old_r2);
                    delta_energy += de;

                }
            }
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
                    if(l == index)
                    {
                        continue;
                    }
                    double new_r2 = min_d2(new_particle.position,S.P[l].position,S.BoxLength);
                    double de = E.morse_potential(new_r2)+E.Debye_Huckel(new_r2);
                    delta_energy += de;

                }
            }
            if(Glauber(delta_energy,gsl_rng_uniform(S.gsl_r)))
            {
                accept += 1.0;
                energy += delta_energy;    
                // Remove from old grid first
                S.G[oldgID].n -= 1;
                S.G[oldgID].plist.remove(index);
                // Add to new grid
                S.G[new_particle.gID].n += 1;
                S.G[new_particle.gID].plist.push_back(index);
                S.P[index] = new_particle;
    
            }
            
        }
        Find_Neighbors();
        Cluster_Particles();
    }
    return accept/double(S.NMOL);
}
//Returns acceptance fraction
// delta is the energy difference
// rand is the random number
bool MC::Glauber(double delta, double rand)
{
    // If the energy difference is greater than 9000 (huge energy difference), the move is not accepted
    if (delta>9000)
        return false;
    else
    {
        // probability of accepting the move is 1/(exp(delta)+1)
        if(1.0/(exp(delta)+1.0)>rand)
            return true;
        else
            return false;
    }
    
}

void MC::Find_Neighbors()
{
    // Clear existing neighbor lists for all particles
    for(int i=0; i<S.NMOL; i++)
    {
        S.P[i].nbr_particles.clear();
    }
    
    // Loop through all grid cells
    for(int gID=0; gID<S.G.size(); gID++)
    {
        // Skip empty grid cells
        if(S.G[gID].plist.empty())
        {
            continue;
        }
        
        // Loop through all particles in this grid cell
        for(int pid : S.G[gID].plist)
        {
            // Loop through all neighboring grid cells
            for(int k=0; k<S.G[gID].nbr.size(); k++)
            {
                int ngID = S.G[gID].nbr[k];
                
                // Skip empty neighboring grid cells
                if(S.G[ngID].plist.empty())
                {
                    continue;
                }
                
                // Loop through all particles in the neighboring grid cell
                for(int npid : S.G[ngID].plist)
                {
                    // Skip self-comparison
                    if(pid == npid || S.P[pid].nbr_particles.find(npid) != S.P[pid].nbr_particles.end())
                    {
                        continue;
                    }
                    
                    // Calculate squared distance between particles
                    double r2 = min_d2(S.P[pid].position, S.P[npid].position, S.BoxLength);
                    
                    // If distance is less than search distance, add to neighbor list
                    if(r2 < S.search2_cm)
                    {
                        S.P[pid].nbr_particles.insert(npid);
                        S.P[npid].nbr_particles.insert(pid);
                    }
                }
            }
        }
    }
}



void MC::Cluster_Particles()
{
    // Clear existing aggregates
    S.Ag.clear();
    
    // Create a vector to track which particles have been processed
    vector<bool> processed(S.NMOL, false);
    
    // Process each particle
    for (int i = 0; i < S.NMOL; i++) {
        // Skip if this particle has already been processed
        if (processed[i]) {
            continue;
        }
        
        // Create a new aggregate
        Aggregate new_aggregate;
        new_aggregate.n = 0;
        new_aggregate.cm = XYZ(0, 0, 0);
        
        // Use breadth-first search to find all connected particles
        queue<int> to_process;
        to_process.push(i);
        processed[i] = true;
        
        while (!to_process.empty()) {
            int current = to_process.front();
            to_process.pop();
            
            // Add current particle to the aggregate
            new_aggregate.plist.push_back(current);
            new_aggregate.n++;
            new_aggregate.cm = new_aggregate.cm + S.P[current].position;
            
            // Update the particle's aggregate ID
            S.P[current].A_ID = S.Ag.size();
            
            // Process all neighbors
            for (int neighbor : S.P[current].nbr_particles) {
                if (!processed[neighbor]) {
                    to_process.push(neighbor);
                    processed[neighbor] = true;
                }
            }
        }
        
        // Calculate the center of mass
        if (new_aggregate.n > 0) {
            new_aggregate.cm.x /= new_aggregate.n;
            new_aggregate.cm.y /= new_aggregate.n;
            new_aggregate.cm.z /= new_aggregate.n;
        }
        
        // Add the new aggregate to the system
        S.Ag.push_back(new_aggregate);
    }

}