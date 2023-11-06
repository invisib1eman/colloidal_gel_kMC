//NANOROD: system.cpp System Class Function Definitions (Revision Date: Oct 27, 2023)
#include "system.h"
void System::Create()
{
    cout<<"Creating System"<<endl;
    
    bool flag;
    
    
    
    int count=0;
  
    for(int i=0; i<NMOL; i++)
    {
        Molecule m;
        m.MOL_ID=i+1;
        m.centre=XYZ(gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L,gsl_rng_uniform(gsl_r)*L-0.5*L);
        m.orientation=RandomRotate(angle_to_quarternion(0,0,0), 1,gsl_rng_uniform(gsl_r),gsl_rng_uniform(gsl_r),gsl_rng_uniform(gsl_r));
        m.UpdateVertices();
        M.push_back(m);
    }
}

void System::WriteMol2(int timestep)
{
    ofstream out;
    out.open("conf.mol2",ios::app);
    out<<timestep<<endl;
    int n_ver=M[0].N_VER;
    
    out<<"@<TRIPOS>MOLECULE"<<endl;
    out<<"Nanorod"<<endl;
    out<<NMOL*(n_ver+1)<<"\t"<<NMOL*n_ver<<endl;
    out<<"SMALL"<<endl;
    out<<"NO_CHARGES"<<endl;

    out<<"@<TRIPOS>ATOM"<<endl;

    string name,type;
    
    int count=0;
    
    XYZ im,shift;
    
    for(int i=0; i<NMOL; i++)
    {
        shift=image(M[i].centre,L)-M[i].centre;
        im=M[i].centre+shift;
        out<<setw(6)<<++count<<"\t"<<"1"<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<"\t"<<"C.3"<<endl;
        
        for(int j=0; j<n_ver; j++)
        {
            im=M[i].ver[j]+shift;
            out<<setw(6)<<++count<<"\t"<<"2"<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<"\t"<<"N.3"<<endl;
        }
    }
          
    out<<"@<TRIPOS>BOND"<<endl;
    
    count=0;
    for(int i=0; i<NMOL; i++)
    {
        for(int j=0; j<n_ver; j++)
        {
            out<<setw(8)<<++count<<"\t"<<setw(8)<<(n_ver+1)*i+1<<"\t"<<setw(8)<<(n_ver+1)*i+j+2<<"\t"<<setw(2)<<"1"<<endl;
        }
    }


    out.close();
    //   cout<<"Mol2 input written in\t"<<FileName<<endl;
    return;
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
    int NT=NMOL*7;
    out<<NT<<endl;
    out<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
    out<<-0.5*L<<"\t"<<0.5*L<<endl;
    out<<-0.5*L<<"\t"<<0.5*L<<endl;
    out<<-0.5*L<<"\t"<<0.5*L<<endl;
    out<<"ITEM: ATOMS index type x y z"<<endl;
    /*int NT=NMOL*7;
    out<<NT<<endl;
    out<<"Time="<<timestep<<endl;*/
    XYZ im;
    XYZ im_centre;
    for(int i=0;i<NMOL;i++)
    {   
        im_centre=image(M[i].centre,L);
        out<<setw(6)<<i*7+1<<"\t"<<1<<"\t"<<setw(8)<<im_centre.x<<"\t"<<setw(8)<<im_centre.y<<"\t"<<setw(8)<<im_centre.z<<endl;
        for(int j=0;j<6;j++)
        {
            im=image(M[i].ver[j],L);
            out<<setw(6)<<i*7+j+2<<"\t"<<2<<"\t"<<setw(8)<<im.x<<"\t"<<setw(8)<<im.y<<"\t"<<setw(8)<<im.z<<endl;
        }
    }
    out.close();
}
