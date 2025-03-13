//NANOROD: system.cpp System Class Function Definitions (Revision Date: Oct 27, 2023)
#include "system.h"

void System::Create()
{
    cout<<NMOL<<endl;
    cout<<"Creating System"<<endl;
    int i,j,k,l,n;
    //Allocate Grid
    Grid g;
    for(k=0; k<NGRID; k++)
    {
        for(j=0; j<NGRID; j++)
        {
          for(i=0; i<NGRID; i++)
          {
              g.g_index=GridIndex_index(i,j,k,NGRID);
              g.nbr.clear();
              g.plist.clear();
              g.cm.x=(double(i)+0.5)*GRIDL-0.5*BoxLength;
              g.cm.y=(double(j)+0.5)*GRIDL-0.5*BoxLength;
              g.cm.z=(double(k)+0.5)*GRIDL-0.5*BoxLength;
              for (int inei=i-1;inei<i+2;inei++)
              {
                for (int jnei=j-1;jnei<j+2;jnei++)
                {
                    for (int knei=k-1;knei<k+2;knei++)
                    {
                        int indexn=GridIndex_index(inei,jnei,knei,NGRID);
                        g.nbr.push_back(indexn);
                        
                    }
                }
              }
              g.n=0;//We fill this later now
              G.push_back(g);
             
          }
        }
    }
   
    //Fill particles Use random permutation
    bool flag = false;
    for(int i=0; i<NMOL; i++)
    {
        Particle p;
        p.P_ID=i;
        // Generate a random position for the particle, make sure there is no overlap with other particles
        flag = true;
        while(flag)
        {
        p.position=XYZ(gsl_rng_uniform(gsl_r)*BoxLength-0.5*BoxLength,gsl_rng_uniform(gsl_r)*BoxLength-0.5*BoxLength,gsl_rng_uniform(gsl_r)*BoxLength-0.5*BoxLength);
        flag = false;
        for(int j=0; j<i; j++)
        {
            if(min_d2(p.position,P[j].position,BoxLength)<cm_L)
            {
                flag = true;
                break;   
            } 
        }
        
        }
        // Calculate the grid index of the particle, and add the particle to the grid
        p.gID=GridIndex_xyz(p.position,NGRID,GRIDL,BoxLength);
        G[p.gID].n+=1;
        G[p.gID].plist.push_back(p.P_ID);
        // Create an aggregate with the particle (each particle is an aggregate at the beginning)
        Aggregate a;
        a.n=1;
        a.plist.push_back(p.P_ID);
        a.rg=1.0;
        a.cm=p.position;
        Ag.push_back(a);
        p.A_ID = i;
        P.push_back(p);
    }
    cout<<"Number of aggregates: "<<Ag.size()<<endl;
    cout<<"Number of particles: "<<P.size()<<endl;
    cout<<"Number of grid: "<<G.size()<<endl;
    int total_particles=0;

    for(int i=0;i<NGRID3;i++)
    {
        /*if(S.G[i].nbr.size()!=nbr_g)
        {
            cout<<"Neighbor Number wrong"<<endl;
            exit(0);
        }*/
        if(G[i].n!=G[i].plist.size())
        {
            cout<<"Init: Grid Number wrong"<<endl;
            exit(0);
        }
        total_particles+=G[i].n;
    }
    if(total_particles!=NMOL)
    {
            cout<<"Init: Total number wrong"<<endl;
            exit(0);
    }
    CreateDump();
}

// void System::WriteMol2(int timestep)
// {
//     ofstream out;
//     out.open("conf.mol2",ios::app);
//     out<<timestep<<endl;
//     int n_ver=M[0].N_VER;
    
//     out<<"@<TRIPOS>MOLECULE"<<endl;
//     out<<"Nanorod"<<endl;
//     out<<NMOL*(n_ver+1)<<"\t"<<NMOL*n_ver+H.size()<<endl;
//     out<<"SMALL"<<endl;
//     out<<"NO_CHARGES"<<endl;

//     out<<"@<TRIPOS>ATOM"<<endl;

//     string name,type;
    
//     int count=0;
    
//     XYZ im,shift;
    
//     for(int i=0; i<NMOL; i++)
//     {
//         shift=image(M[i].centre,L)-M[i].centre;
//         im=M[i].centre+shift;
//         out<<setw(6)<<++count<<"\t"<<"1"<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<"\t"<<"C.3"<<endl;
        
//         for(int j=0; j<n_ver; j++)
//         {
//             im=M[i].ver[j]+shift;
//             out<<setw(6)<<++count<<"\t"<<"2"<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<"\t"<<"N.3"<<endl;
//         }
//     }
          
//     out<<"@<TRIPOS>BOND"<<endl;
    
//     count=0;
//     for(int i=0; i<NMOL; i++)
//     {
//         for(int j=0; j<n_ver; j++)
//         {
//             out<<setw(8)<<++count<<"\t"<<setw(8)<<(n_ver+1)*i+1<<"\t"<<setw(8)<<(n_ver+1)*i+j+2<<"\t"<<setw(2)<<"1"<<endl;
//         }
//     }
//     list<hbond>::iterator it;
//     for(it=H.begin();it!=H.end();it++)
//     {
//         hbond writebond=*it;
//         out<<setw(8)<<++count<<"\t"<<setw(8)<<(n_ver+1)*writebond.M1+writebond.arm1+1<<"\t"<<setw(8)<<(n_ver+1)*writebond.M2+writebond.arm2+1<<"\t"<<setw(2)<<"1"<<endl;
            
        
        
//     }
        
    

//     out.close();
//     //   cout<<"Mol2 input written in\t"<<FileName<<endl;
//     return;
// }
void System::CreateDump()
{
    string FileName = dump_file_name;
    ofstream out;
    out.open(FileName,ios::trunc);
    out.close();
}
void System::WriteDump(int timestep)
{
    string FileName = dump_file_name;
    ofstream out;
    out.open(FileName,ios::app);
    
    out<<"ITEM: TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<"ITEM: NUMBER OF ATOMS"<<endl;
    int NT=NMOL*1;
    out<<NT<<endl;
    out<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
    out<<-0.5*BoxLength<<"\t"<<0.5*BoxLength<<endl;
    out<<-0.5*BoxLength<<"\t"<<0.5*BoxLength<<endl;
    out<<-0.5*BoxLength<<"\t"<<0.5*BoxLength<<endl;
    out<<"ITEM: ATOMS index type x y z"<<endl;
    /*int NT=NMOL*7;
    out<<NT<<endl;
    out<<"Time="<<timestep<<endl;*/
    XYZ im;
    XYZ im_pos;
    for(int i=0;i<NMOL;i++)
    {   
        im_pos=image(P[i].position,BoxLength);
        out<<setw(6)<<i+1<<"\t"<<1<<"\t"<<setw(8)<<im_pos.x<<"\t"<<setw(8)<<im_pos.y<<"\t"<<setw(8)<<im_pos.z<<endl;
        
    }
    out.close();
}
void System::WriteData(int timestep)
{
    string FileName = data_file_name;
    ofstream out;
    out.open(FileName,ios::trunc);
    
    // Header section
    out << "LAMMPS data file from NANOROD simulation at timestep " << timestep << endl;
    out << endl;
    out << NMOL << " atoms" << endl;
    out << "1 atom types" << endl;
    out << endl;
    
    // Box dimensions
    out << -0.5*BoxLength << " " << 0.5*BoxLength << " xlo xhi" << endl;
    out << -0.5*BoxLength << " " << 0.5*BoxLength << " ylo yhi" << endl; 
    out << -0.5*BoxLength << " " << 0.5*BoxLength << " zlo zhi" << endl;
    out << endl;

    // Atoms section
    out << "Atoms" << endl;
    out << endl;

    XYZ im_centre;
    for(int i=0; i<NMOL; i++)
    {   
        im_centre = image(P[i].position,BoxLength);
        // Format: atom-ID atom-type x y z
        out << i+1 << " " << i+1 << " " << 1 << " " 
            << fixed << setprecision(6) << im_centre.x << " " 
            << im_centre.y << " " 
            << im_centre.z << endl;
    }
    out.close();
}
// void System::WriteBond(int timestep)
// {
//     ofstream out;
//     out.open("Bondlist.txt",ios::app);
//     out<<"TIMESTEP"<<endl;
//     out<<timestep<<endl;
//     out<<setw(12)<<"molecule1"<<"\t"<<setw(12)<<"molecule2"<<"\t"<<setw(12)<<"arm1"<<"\t"<<setw(8)<<"arm2"<<endl;
//     list<hbond>::iterator it;
//     for(it=H.begin();it!=H.end();it++)
//     {
//         hbond writebond=*it;
//         out<<setw(12)<<writebond.M1<<setw(12)<<writebond.M2<<setw(12)<<writebond.arm1<<setw(12)<<writebond.arm2<<endl;
            
        
        
//     }
//     out.close();
// }
// void System::WriteOrientation(int timestep)
// {
//     ofstream out;
//     out.open("Orientation.txt",ios::app);
//     out<<"TIMESTEP"<<endl;
//     out<<timestep<<endl;
//     out<<setw(12)<<"molecule1"<<"\t"<<setw(12)<<"molecule2"<<"\t"<<setw(12)<<"arm1"<<"\t"<<setw(8)<<"arm2"<<endl;
//     for(int j=0;j<NMOL;j++)
//     {
        
//         out<<setw(12)<<M[j].orientation.w<<setw(12)<<M[j].orientation.x<<setw(12)<<M[j].orientation.y<<setw(12)<<M[j].orientation.z<<endl;
            
        
    
//     }
//     out.close();

// }
// void System::WriteGrid(int timestep)
// {
//     ofstream out;
//     out.open("Grid.txt",ios::app);
//     out<<"TIMESTEP"<<endl;
//     out<<timestep<<endl;
//     for(int i=0;i<NGRID3;i++)
//     {
//      out<<setw(12)<<G[i].g_index<<setw(12)<<G[i].cm.x<<setw(12)<<G[i].cm.y<<setw(12)<<G[i].cm.z<<setw(12)<<G[i].n<<setw(8)<<endl;
//               out<<"neighborlist"<<endl;
//               for(int l=0;l<G[i].nbr.size();l++)
//               {
//                 out<<G[i].nbr[l]<<endl;
//               }
              
//     }
//     out.close();
// }
void System::UpdateGrid()
{
    
    int i,j,k,l,n;
    //Allocate Grid
    int g_index;
    for(k=0; k<NGRID; k++)
    {
        for(j=0; j<NGRID; j++)
        {
          for(i=0; i<NGRID; i++)
          {
              
              g_index=GridIndex_index(i,j,k,NGRID);
              
             
              G[g_index].cm.x=(double(i)+0.5)*GRIDL-0.5*BoxLength;
              G[g_index].cm.y=(double(j)+0.5)*GRIDL-0.5*BoxLength;
              G[g_index].cm.z=(double(k)+0.5)*GRIDL-0.5*BoxLength;
              for (int inei=0;inei<3;inei++)
              {
                for (int jnei=0;jnei<3;jnei++)
                {
                    for (int knei=0;knei<3;knei++)
                    {
                        int indexn=GridIndex_index(inei+i-1,jnei+j-1,knei+k-1,NGRID);
                        G[g_index].nbr[inei*9+jnei*3+knei]=indexn;
                        
                    }
                }
              }
              
             
          }
        }
    }
}
void System::WriteRestart()
{
    string filename = restart_file_name;
    ofstream out(filename);
    if (!out.is_open()) {
        cout << "Error: Could not open restart file" << endl;
        exit(1);
    }

    // Write number of particles
    out << NMOL << endl;

    // Write box length
    out << BoxLength << endl;

    // Write particle positions and properties
    for (int i = 0; i < NMOL; i++) {
        out << setw(12) << P[i].position.x 
            << setw(12) << P[i].position.y
            << setw(12) << P[i].position.z
            << setw(12) << P[i].P_ID 
            << setw(12) << P[i].A_ID
            << setw(12) << P[i].gID << endl;
    }

    // Write grid information
    out << NGRID << endl;
    out << NGRID3 << endl;
    out << GRIDL << endl;

    // Write grid cell data
    for (int i = 0; i < NGRID3; i++) {
        out << setw(12) << G[i].cm.x
            << setw(12) << G[i].cm.y 
            << setw(12) << G[i].cm.z;
        
        // Write particle list for each grid
        out << setw(12) << G[i].plist.size();
        for (auto p : G[i].plist) {
            out << setw(12) << p;
        }
        out << endl;
    }

    // Write aggregate information
    out << Ag.size() << endl;
    for (int i = 0; i < Ag.size(); i++) {
        out << setw(12) << Ag[i].n  // Number of particles in aggregate
            << setw(12) << Ag[i].cm.x
            << setw(12) << Ag[i].cm.y
            << setw(12) << Ag[i].cm.z << endl;
        
        // Write particle IDs in this aggregate
        for (int j = 0; j < Ag[i].n; j++) {
            out << setw(12) << Ag[i].plist[j];
        }
        out << endl;
    }

    out.close();
}
void System::ReadRestart()
{
    string filename = read_restart_file_name;
    ifstream in(filename);
    if (!in.is_open()) {
        cout << "Error: Could not open restart file" << endl;
        exit(1);
    }

    // Read system parameters
    in >> NMOL;
    P.resize(NMOL);

    // Read particle data
    for (int i = 0; i < NMOL; i++) {
        in >> P[i].position.x 
           >> P[i].position.y
           >> P[i].position.z
           >> P[i].P_ID
           >> P[i].A_ID
           >> P[i].gID;
    }

    // Read grid parameters
    in >> NGRID;
    in >> NGRID3;
    in >> GRIDL;

    // Read grid cell data
    G.resize(NGRID3);
    for (int i = 0; i < NGRID3; i++) {
        in >> G[i].cm.x
           >> G[i].cm.y
           >> G[i].cm.z;
        
        int plist_size;
        in >> plist_size;
        
        // Read particle list for each grid
        for (int j = 0; j < plist_size; j++) {
            int p;
            in >> p;
            G[i].plist.push_back(p);
        }
    }

    // Read aggregate data
    int n_aggregates;
    in >> n_aggregates;
    Ag.resize(n_aggregates);
    
    for (int i = 0; i < n_aggregates; i++) {
        in >> Ag[i].n
           >> Ag[i].cm.x
           >> Ag[i].cm.y
           >> Ag[i].cm.z;

        // Read particle IDs in this aggregate
        Ag[i].plist.resize(Ag[i].n);
        for (int j = 0; j < Ag[i].n; j++) {
            in >> Ag[i].plist[j];
        }
    }

    in.close();
    
}


