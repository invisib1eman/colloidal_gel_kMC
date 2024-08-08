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
              g.cm.x=(double(i)+0.5)*GRIDL-0.5*L;
              g.cm.y=(double(j)+0.5)*GRIDL-0.5*L;
              g.cm.z=(double(k)+0.5)*GRIDL-0.5*L;
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
    bool flag;
    
    
    int count=0;
    /*Molecule m1;//initial test for 2 particles
    m1.MOL_ID=0;
    m1.centre=XYZ(0,0,0);
    m1.gID=GridIndex_xyz(m1.centre,NGRID,GRIDL,L);
    G[m1.gID].n+=1;
    G[m1.gID].plist.push_back(m1.MOL_ID);
    m1.UpdateVertices();
    m1.nbonds=1;
    m1.hbond_list.push_back(hbond(0,1,0,4));
    m1.vertype[0]='I';
    M.push_back(m1);
    Molecule m2;
    m2.MOL_ID=1;
    m2.centre=XYZ(2.26,0,1.06);
    m2.gID=GridIndex_xyz(m2.centre,NGRID,GRIDL,L);
    G[m2.gID].n+=1;
    G[m2.gID].plist.push_back(m2.MOL_ID);
    m2.orientation=angle_to_quarternion(M_PI/3,1,0);
    m2.UpdateVertices();
    m2.nbonds=1;
    m2.hbond_list.push_back(hbond(1,0,4,0));
    m2.vertype[4]='I';
    M.push_back(m2);
    NMOL=2;*/
    for(int i=0; i<NMOL; i++)
    {
        Particle p;
        p.P_ID=i;
        p.position=XYZ(gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L);
        p.gID=GridIndex_xyz(p.position,NGRID,GRIDL,L);
        G[p.gID].n+=1;
        G[p.gID].plist.push_back(p.P_ID);
        P.push_back(p);
    }
    cout<<"initialized"<<endl;

}



void System::WriteDump(int timestep)
{
    char FileName[200];
    sprintf(FileName,"%s_Dump.lammpstrj",Description.c_str());
    ofstream out;
    out.open(FileName,ios::app);
    
    out<<"ITEM: TIMESTEP"<<endl;
    out<<timestep<<endl;
    out<<"ITEM: NUMBER OF ATOMS"<<endl;
    out<<NMOL<<endl;
    out<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
    out<<-0.5*L<<"\t"<<0.5*L<<endl;
    out<<-0.5*L<<"\t"<<0.5*L<<endl;
    out<<-0.5*L<<"\t"<<0.5*L<<endl;
    out<<"ITEM: ATOMS index type x y z"<<endl;
    /*int NT=NMOL*7;
    out<<NT<<endl;
    out<<"Time="<<timestep<<endl;*/
    XYZ im;
    XYZ im_pos;
    for(int i=0;i<NMOL;i++)
    {   
        im_pos=image(P[i].position,L);
        out<<setw(6)<<i+1<<"\t"<<1<<"\t"<<setw(8)<<im_pos.x<<"\t"<<setw(8)<<im_pos.y<<"\t"<<setw(8)<<im_pos.z<<endl;
        
    }
    out.close();
}


void System::WriteGrid(int timestep)
{
    ofstream out;
    out.open("Grid.txt",ios::app);
    out<<"TIMESTEP"<<endl;
    out<<timestep<<endl;
    for(int i=0;i<NGRID3;i++)
    {
     out<<setw(12)<<G[i].g_index<<setw(12)<<G[i].cm.x<<setw(12)<<G[i].cm.y<<setw(12)<<G[i].cm.z<<setw(12)<<G[i].n<<setw(8)<<endl;
              out<<"neighborlist"<<endl;
              for(int l=0;l<G[i].nbr.size();l++)
              {
                out<<G[i].nbr[l]<<endl;
              }
              
    }
    out.close();
}
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
              
             
              G[g_index].cm.x=(double(i)+0.5)*GRIDL-0.5*L;
              G[g_index].cm.y=(double(j)+0.5)*GRIDL-0.5*L;
              G[g_index].cm.z=(double(k)+0.5)*GRIDL-0.5*L;
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


