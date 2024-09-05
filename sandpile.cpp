/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <ttp://www.gnu.org/licenses/>  *
 ************************************************************************/

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;

struct UserData
{

    
    double height;
    double R;
    double rho_s;
    double acc;
    double tp;
    double ADDS=0.0;
    double Kn;
    double Kt;
    double Gn;
    double MUP;
    std::ofstream      oss_ss;       ///< file for stress strain data
};


void Report (DEM::Domain & d2, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (d2.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",d2.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "Nc" << Util::_8s << "Nsc \n";
    }
    if (!d2.Finished) 
    {
        size_t Nc = 0;
        size_t Nsc = 0;
        for (size_t i=0; i<d2.CInteractons.Size(); i++)
        {
            Nc += d2.CInteractons[i]->Nc;    //Number of contact
            Nsc += d2.CInteractons[i]->Nsc;  // Number of sliding contact
        }
        dat.oss_ss << Util::_10_6 << d2.Time << Util::_8s << Nc << Util::_8s << Nsc << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

void Setup (DEM::Domain & d2, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double tp=pow(2.0*2.0*dat.R/abs(dat.acc), 0.5);
    double Nstep=ceil(tp*dat.tp/d2.Dt);
    //cout <<"datacc"<<dat.acc<<endl;
    //cout <<"datR"<<dat.R<<endl;
    //cout <<"Nstep"<<Nstep<<endl;
    //cout <<"Tp"<<tp<<endl;
    //Nstep=15;
    if ( fmod(d2.iter,  Nstep)-0.0 ==0.0)
    {
        dat.ADDS=dat.ADDS+1;
        cout <<"Add Sphere"<<dat.ADDS<<endl;
        cout <<"Nstep"<<Nstep<<endl;
        cout <<"Tp"<<tp<<endl;
        cout <<"fmod(d2.iter,  Nstep)"<<fmod(d2.iter,  Nstep)<<endl;
        
        d2.AddSphere(-1,Vec3_t(0.0,0.0,dat.height),dat.R,dat.rho_s);
        d2.Particles[d2.Particles.Size()-1]->Props.Kn = dat.Kn; // normal stiffness
        d2.Particles[d2.Particles.Size()-1]->Props.Kt = dat.Kt; // trangential stiffness
        d2.Particles[d2.Particles.Size()-1]->Props.Gn = dat.Gn; // restitution coefficient
        d2.Particles[d2.Particles.Size()-1]->Props.Mu = dat.MUP; // frictional coefficient
        d2.Particles[d2.Particles.Size()-1]->Ff = d2.Particles[d2.Particles.Size()-1]->Props.m*Vec3_t(0.0,0.0,dat.acc);
        d2.Particles[d2.Particles.Size()-1]->v = Vec3_t(0.0,0.0,0.0);
        d2.Particles[d2.Particles.Size()-1]->InitializeVelocity(d2.Dt);
        cout <<"Add V"<<d2.Particles[d2.Particles.Size()-1]->v<<endl;
    }
    #pragma omp parallel for schedule(static) num_threads(d2.Nproc)
    for (size_t np=0;np<d2.Particles.Size();np++)
    {
        d2.Particles[np]->InitializeVelocity(d2.Dt);
    }
    d2.UpdateContacts();
    //cout <<"v"<< vmod<<endl;
 
}



int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////

    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    //number of threads
    size_t Nproc = 1; 
    if (argc==3) Nproc=atoi(argv[2]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());
    
    String ptype;
    bool   Render = true;
    double L_x       ;
    double L_y       ;
    double L_z       ;
    double nu        ;
    double dx        ;
    double dt        ;
    double Dp        ;
    double R         ;   
    double Tf1       ;
    double Tf2       ;
    double dtOut1    ;
    double dtOut2    ;
    double height    ;
    double rho_s     ;
    double rho_l     ;
    double Mu        ;
    double Muw       ;
    double Beta      ;        //
	double Kn        ;         // Normal stiffness
	double Kt        ;         // Tangential stiffness
	double Gn        ;         // Normal dissipative coefficient
	double ratiots   ;     // Coefficient to modify the time step
	double ratiokn   ;
    double R_random  ;
    double Gnw       ;
    double sizedis    ;
    double MUP       ;
    double Ang       ;
    double fraction  ;
    int    Nump      ;
    double Amp       ;
    double tp      ;



    
    {
        infile >> ptype;           infile.ignore(200,'\n');
        infile >> Render;          infile.ignore(200,'\n');
        infile >> L_x;             infile.ignore(200,'\n');
        infile >> L_y;             infile.ignore(200,'\n');
        infile >> L_z;             infile.ignore(200,'\n');
        infile >> nu;              infile.ignore(200,'\n');
        infile >> dx;              infile.ignore(200,'\n');
        infile >> dt;              infile.ignore(200,'\n');
        infile >> Dp;              infile.ignore(200,'\n');
        infile >> R;               infile.ignore(200,'\n');
        infile >> Tf1;             infile.ignore(200,'\n');
        infile >> Tf2;             infile.ignore(200,'\n');
        infile >> dtOut1;          infile.ignore(200,'\n');
        infile >> dtOut2;          infile.ignore(200,'\n');
        infile >> height;          infile.ignore(200,'\n');
        infile >> rho_s;           infile.ignore(200,'\n');
        infile >> rho_l;           infile.ignore(200,'\n');
        infile >> Mu;              infile.ignore(200,'\n');
        infile >> Muw;             infile.ignore(200,'\n');
        infile >> Beta;            infile.ignore(200,'\n');
        infile >> Kn;              infile.ignore(200,'\n');
        infile >> Gn;              infile.ignore(200,'\n');
        infile >> ratiots;         infile.ignore(200,'\n');
        infile >> ratiokn;         infile.ignore(200,'\n');
        infile >> R_random;        infile.ignore(200,'\n');
        infile >> Gnw;             infile.ignore(200,'\n');
        infile >> sizedis;          infile.ignore(200,'\n');
        infile >> MUP;             infile.ignore(200,'\n');
        infile >> Ang;             infile.ignore(200,'\n');
        infile >> fraction;             infile.ignore(200,'\n');
        infile >> Nump;             infile.ignore(200,'\n');
        infile >> Amp;             infile.ignore(200,'\n');
        infile >> tp;             infile.ignore(200,'\n');

    }

    Kt=ratiokn*Kn;
  

    double acc=Dp * (rho_s - rho_l) / rho_s;
    double ang=Ang*M_PI/180.0;
    cout <<"tp"<< tp<<endl;
    cout <<"L_z"<< L_z<<endl;
    double Sphere_size;
    double dtdem;
    Vec3_t Xmin0;
    Vec3_t Xmax0;
    Vec3_t Center;
    double R_base=9.2;
    Array<int> delpar0;
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    double thichness=0.1;
    


    DEM::Domain d2;
    
    UserData dat;
    d2.UserData = &dat;
    dat.R=R;
    dat.rho_s=rho_s;
    dat.acc=acc;
    dat.height=height;
    dat.tp=tp;
    dat.MUP=MUP;
    dat.Kn=Kn;
    dat.Kt=Kt;
    Center = Vec3_t(0.0,0.0,0.0);

    Vec3_t Xmin0base=Vec3_t(Center(0)-R_base, Center(1)-R_base    , 0.0);
    Vec3_t Xmax0base=Vec3_t(Center(0)+R_base, Center(1)+R_base    , R);
    
    d2.AddSphere(-1,Vec3_t(0.0,0.0,height-16.0*R),R,rho_s);
    double scaleL=1.2;
    d2.GenSpheresBox (-15, Xmin0base, Xmax0base, R/2.0, rho_s, "HCP",  1234, 1.0, sizedis);
    d2.AddPlane (-16, Vec3_t(Center(0),Center(1),Center(2)-L_z), R, scaleL*R_base*2.0, scaleL*R_base*2.0, 1.0, 0.0, &axis1);
    d2.AddPlane (-17, Vec3_t(Center(0)-R_base*scaleL,Center(1),Center(2)), R, L_z+height, scaleL*R_base*2.0, 1.0, ang, &axis1);
    d2.AddPlane (-18, Vec3_t(Center(0)+R_base*scaleL,Center(1),Center(2)), R, L_z+height, scaleL*R_base*2.0, 1.0, ang, &axis1);
    d2.AddPlane (-19, Vec3_t(Center(0),Center(1)-R_base*scaleL,Center(2)), R, scaleL*R_base*2.0, L_z+height, 1.0, M_PI/2.0, &axis0);
    d2.AddPlane (-20, Vec3_t(Center(0),Center(1)+R_base*scaleL,Center(2)), R, scaleL*R_base*2.0, L_z+height, 1.0, M_PI/2.0, &axis0);
    d2.GetParticle(-16)->FixVeloc(); 
    d2.GetParticle(-17)->FixVeloc(); 
    d2.GetParticle(-18)->FixVeloc(); 
    d2.GetParticle(-19)->FixVeloc(); 
    d2.GetParticle(-20)->FixVeloc(); 
    size_t countdel = 0;

    for (size_t np=0;np<d2.Particles.Size();np++)
    {

        if (d2.Particles[np]->Tag == -15&& pow((d2.Particles[np]->x(0)-Center(0)),2)+pow((d2.Particles[np]->x(1)-Center(1)),2)+pow((d2.Particles[np]->x(2)-Center(2)),2)>=R_base*R_base)
        {
            countdel = countdel+1;
            d2.Particles[np]->Tag = 10;
        }

    }

    Array<int> delpar1;
    delpar1.Push(10);
    d2.DelParticles(delpar1);


    for (size_t np=0;np<d2.Particles.Size();np++)
    {
        if (d2.Particles[np]->Tag == -1)
        {
        
            d2.Particles[np]->Ff = d2.Particles[np]->Props.m*Vec3_t(0.0,0.0,acc);
            d2.Particles[np]->Props.Kn = Kn; // normal stiffness
            d2.Particles[np]->Props.Kt = Kt; // trangential stiffness
            d2.Particles[np]->Props.Gn = Gn; // restitution coefficient
            d2.Particles[np]->Props.Mu = MUP; // frictional coefficient
        }

        if (d2.Particles[np]->Tag == -15)
        {
            d2.Particles[np]->v==Vec3_t(0.0,0.0,0.0);
            d2.Particles[np]->FixVeloc();
            d2.Particles[np]->Props.Kn = Kn; // normal stiffness
            d2.Particles[np]->Props.Kt = Kt; // trangential stiffness
            d2.Particles[np]->Props.Gn = Gn; // restitution coefficient
            d2.Particles[np]->Props.Mu = Muw; // frictional coefficient
        }
        if (d2.Particles[np]->Tag < -15)
        {
            d2.Particles[np]->v==Vec3_t(0.0,0.0,0.0);
      
            d2.Particles[np]->Props.Kn = Kn; // normal stiffness
            d2.Particles[np]->Props.Kt = Kt; // trangential stiffness
            d2.Particles[np]->Props.Gn = Gn; // restitution coefficient
            d2.Particles[np]->Props.Mu = Muw; // frictional coefficient
        }
    }

    dtdem = ratiots*d2.CriticalDt(); //Calculating time step
    cout <<"dtdem"<< dtdem<<endl;
    cout <<"ratiots"<< ratiots<<endl;
    cout <<"CACU"<< d2.CriticalDt()<<endl;

    d2.Alpha = R; //Verlet distance
    d2.Solve(/*tf*/Tf2, dtdem, /*dtOut*/dtOut2, Setup, Report, "sandpile", 2, Nproc);
    d2.Save("Stage_cpu");

}
MECHSYS_CATCH
