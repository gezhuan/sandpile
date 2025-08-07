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
#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>

using std::cout;
using std::endl;

struct UserData
{

    
    double height;
    double R;
    double rho_s;
    double acc;
    double tp;
    int    ADDS=0;
    double Kn;
    double Kt;
    double Gn;
    double MUP;
    double timeperadd;
    double R_source;
    double w_th;
    String ptype;
    int Nump;
    int Anumper;
    std::ofstream      oss_ss;       ///< file for stress strain data
};

struct Circle {
    double x;
    double y;
    double r;
};

// 检查两个圆是否重叠
bool circlesOverlap(const Circle& c1, const Circle& c2) {
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    double distance = std::sqrt(dx*dx + dy*dy);
    return distance < (c1.r + c2.r);
}

// 检查圆是否在大圆内
bool isInsideBigCircle(const Circle& small, const Circle& big) {
    double dx = small.x - big.x;
    double dy = small.y - big.y;
    double distance = std::sqrt(dx*dx + dy*dy);
    return distance <= (big.r - small.r);
}

// 生成N个不重叠的小圆
std::vector<Circle> generateNonOverlappingCircles(
    double bigR, double bigX, double bigY, double smallR, 
    int N, int maxAttemptsPerCircle = 1000, int maxTotalAttempts = 10000) {
    
    // 准备随机数生成器
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> angleDist(0.0, 2 * M_PI);
    std::uniform_real_distribution<double> radiusDist(0.0, 1.0);
    
    Circle bigCircle = {bigX, bigY, bigR};
    std::vector<Circle> smallCircles;
    int totalAttempts = 0;
    
    for (int i = 0; i < N && totalAttempts < maxTotalAttempts; ) {
        int attempts = 0;
        bool placed = false;
        
        while (!placed && attempts < maxAttemptsPerCircle && totalAttempts < maxTotalAttempts) {
            // 生成随机角度和半径（使用平方根确保均匀分布）
            double angle = angleDist(generator);
            double radius = std::sqrt(radiusDist(generator)) * (bigR - smallR);
            
            // 计算坐标
            Circle newCircle = {
                bigX + radius * std::cos(angle),
                bigY + radius * std::sin(angle),
                smallR
            };
            
            // 检查是否在大圆内且不与任何已有小圆重叠
            bool valid = isInsideBigCircle(newCircle, bigCircle);
            for (const auto& circle : smallCircles) {
                if (circlesOverlap(newCircle, circle)) {
                    valid = false;
                    break;
                }
            }
            
            if (valid) {
                smallCircles.push_back(newCircle);
                placed = true;
                i++;  // 只有成功放置才增加i
            }
            
            attempts++;
            totalAttempts++;
        }
        
        if (!placed) {
            std::cerr << "警告: 无法放置第 " << i+1 << " 个小圆，已尝试 " 
                      << attempts << " 次\n";
            break;  // 如果无法放置当前小圆，直接退出
        }
    }
    
    if (smallCircles.size() < N) {
        std::cerr << "警告: 仅成功放置 " << smallCircles.size() 
                  << " 个中的 " << N << " 个小圆\n";
    }
    
    return smallCircles;
}


void Setup (DEM::Domain & d2, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double tp=pow(2.0*2.0*dat.R/abs(dat.acc), 0.5);
    double Nstep=ceil(tp*dat.tp/d2.Dt);
    dat.timeperadd=tp;
    Vec3_t axis0(OrthoSys::e0); // rotation of face

    //cout <<"datacc"<<dat.acc<<endl;
    //cout <<"datR"<<dat.R<<endl;
    //cout <<"Nstep"<<Nstep<<endl;
    //cout <<"Tp"<<tp<<endl;
    //Nstep=15;
    if ( fmod(d2.iter,  Nstep)-0.0 ==0.0&&dat.ADDS< dat.Nump)
    {
        
        //cout <<"Add Sphere"<<dat.ADDS<<endl;
        //cout <<"Nstep"<<Nstep<<endl;
        //cout <<"dat.Gn"<<dat.Gn<<endl;
        //cout <<"Tp"<<tp<<endl;
        //cout <<"fmod(d2.iter,  Nstep)"<<fmod(d2.iter,  Nstep)<<endl;


        double bigR = dat.R_source*dat.R;    // 大圆半径
        double bigX = 0.0;     // 大圆中心x坐标
        double bigY = 0.0;     // 大圆中心y坐标
        double smallR = dat.R;   // 小圆半径
 
        
        auto circles = generateNonOverlappingCircles(bigR, bigX, bigY, smallR, dat.Anumper);
        for (size_t np=0;np<circles.size();np++)
        {
            dat.ADDS=dat.ADDS+1;
            //d2.AddSphere(-1,Vec3_t(circles[np].x,circles[np].y,dat.height),dat.R,dat.rho_s);
            if (dat.ptype=="sphere") d2.AddSphere (-1,Vec3_t(circles[np].x,circles[np].y,dat.height),dat.R,dat.rho_s);
            if (dat.ptype=="cube")   d2.AddCube   (-1, Vec3_t(circles[np].x,circles[np].y,dat.height),dat.w_th,0.5*dat.R,dat.rho_s);
            if (dat.ptype=="tetra")  d2.AddTetra  (-1, Vec3_t(circles[np].x,circles[np].y,dat.height),dat.w_th,0.5*dat.R,dat.rho_s);
            if (dat.ptype=="rice")  d2.AddRice   (-1, Vec3_t(circles[np].x,circles[np].y,dat.height), dat.w_th, dat.R*2.0, dat.rho_s, M_PI/2.0, &axis0);
        }
        

        for (size_t np=0;np<dat.Anumper;np++)
        {
            d2.Particles[d2.Particles.Size()-np-1]->Props.Kn = dat.Kn; // normal stiffness
            d2.Particles[d2.Particles.Size()-np-1]->Props.Kt = dat.Kt; // trangential stiffness
            d2.Particles[d2.Particles.Size()-np-1]->Props.Gn = dat.Gn; // restitution coefficient
            d2.Particles[d2.Particles.Size()-np-1]->Props.Mu = dat.MUP; // frictional coefficient
            d2.Particles[d2.Particles.Size()-np-1]->Ff = d2.Particles[d2.Particles.Size()-np-1]->Props.m*Vec3_t(0.0,0.0,dat.acc);
            d2.Particles[d2.Particles.Size()-np-1]->v = Vec3_t(0.0,0.0,0.0);
            d2.Particles[d2.Particles.Size()-np-1]->InitializeVelocity(d2.Dt);
            //cout <<"Add V"<<d2.Particles[d2.Particles.Size()-1]->v<<endl;
            d2.FreePar.Push(d2.Particles.Size()-np-1);
        }

        
        d2.UpdateContacts();
    //std::cout << d2.Interactons.Size() << std::endl;
    }

    
    /*#pragma omp parallel for schedule(static) num_threads(d2.Nproc)
    for (size_t np=0;np<d2.Particles.Size();np++)
    {
        if (d2.Particles[np]->Tag == -1&& d2.Particles[np]->x(2)-Center(2)<=3.0*dat.R)
        {
            //countdel = countdel+1;
            DelParticlesIdx(np);
            //d2.Particles[np]->Tag = 10;
        }
    }*/
    //cout <<"v"<< vmod<<endl;
 
}

void Report (DEM::Domain & d2, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (d2.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",d2.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time"<< Util::_8s << "AddN" << Util::_8s << "Tf" << Util::_8s << "Nc" << Util::_8s << "Nsc \n";
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
        dat.oss_ss << Util::_10_6 << d2.Time << Util::_8s << dat.ADDS << Util::_8s << dat.timeperadd<< Util::_8s << Nc << Util::_8s << Nsc << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////

    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    //number of threads
    size_t Nproc = 1; 
    if (argc>=3) Nproc=atoi(argv[2]);
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
    double Tf2       ;
    double dtOut2    ;
    double height    ;
    double rho_s     ;
    double w_th      ;
    double Muw       ;
    double Beta      ;        //
	double Kn        ;         // Normal stiffness
	double Kt        ;         // Tangential stiffness
	double Gn        ;         // Normal dissipative coefficient
    double Gnb       ;         // Normal dissipative coefficient
	double ratiots   ;     // Coefficient to modify the time step
	double ratiokn   ;
    double Rp_base   ;
    double Gnw       ;
    double sizedis   ;
    double MUP       ;
    double Ang       ;
    double fraction  ;
    double R_source  ;
    int    Nump      ;
    double tp        ;
    int    Anumper   ;
    String contactlaw;


    
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
        infile >> Tf2;             infile.ignore(200,'\n');
        infile >> dtOut2;          infile.ignore(200,'\n');
        infile >> height;          infile.ignore(200,'\n');
        infile >> rho_s;           infile.ignore(200,'\n');
        infile >> w_th;             infile.ignore(200,'\n');
        infile >> Muw;             infile.ignore(200,'\n');
        infile >> Beta;            infile.ignore(200,'\n');
        infile >> Kn;              infile.ignore(200,'\n');
        infile >> Gn;              infile.ignore(200,'\n');
        infile >> Gnb;             infile.ignore(200,'\n');
        infile >> ratiots;         infile.ignore(200,'\n');
        infile >> ratiokn;         infile.ignore(200,'\n');
        infile >> Rp_base;         infile.ignore(200,'\n');
        infile >> Gnw;             infile.ignore(200,'\n');
        infile >> sizedis;         infile.ignore(200,'\n');
        infile >> MUP;             infile.ignore(200,'\n');
        infile >> Ang;             infile.ignore(200,'\n');
        infile >> fraction;        infile.ignore(200,'\n');
        infile >> R_source;        infile.ignore(200,'\n');
        infile >> Nump;            infile.ignore(200,'\n');
        infile >> tp;              infile.ignore(200,'\n');
        infile >> Anumper;         infile.ignore(200,'\n');
        infile >> contactlaw;   infile.ignore(200,'\n');

        

    }

    Kt=ratiokn*Kn;
  

    double acc=Dp ;
    double ang=Ang*M_PI/180.0;
    cout <<"tp"<< tp<<endl;
    cout <<"Rp_base"<< Rp_base<<endl;
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
    

    size_t cl = 0;
    if (contactlaw=="hertz") cl=1;
    cout <<"contactlaw"<< contactlaw <<endl;
    UserData dat;
    DEM::Domain d2(&dat,cl);


    d2.UserData = &dat;
    dat.R=R;
    dat.rho_s=rho_s;
    dat.acc=acc;
    dat.height=height;
    dat.tp=tp;
    dat.MUP=MUP;
    dat.Kn=Kn;
    dat.Kt=Kt;
    dat.Gn=Gn;
    dat.Nump=Nump;
    dat.Anumper=Anumper;
    dat.R_source=R_source;
    dat.w_th=w_th;
    dat.ptype=ptype;
    Center = Vec3_t(0.0,0.0,0.0);

    Vec3_t Xmin0base=Vec3_t(Center(0)-1.2*R_base, Center(1)-1.2*R_base    , 0.0);
    Vec3_t Xmax0base=Vec3_t(Center(0)+1.2*R_base, Center(1)+1.2*R_base    , Rp_base*2.0);
 
    if (ptype=="sphere") d2.AddSphere(-1,Vec3_t(0.0,0.0,height-16.0*R),R,rho_s);
    if (ptype=="cube") d2.AddCube (-1, Vec3_t(0.0,0.0,height-16.0*R),w_th,0.5*R,rho_s);
    if (ptype=="tetra") d2.AddTetra (-1, Vec3_t(0.0,0.0,height-16.0*R),w_th,0.5*R,rho_s);
    if (ptype=="rice") d2.AddRice   (-1, Vec3_t(0.0,0.0,height-16.0*R), w_th, R*2.0, rho_s, M_PI/2.0, &axis0);
    double scaleL=2.0;
    d2.GenSpheresBox (-15, Xmin0base, Xmax0base, Rp_base, rho_s, "HCP",  1234, 1.0, sizedis);
    cout <<"Number of Particles"<< d2.Particles.Size()<<endl;
    d2.AddPlane (-16, Vec3_t(Center(0),Center(1),Center(2)-L_z), R, scaleL*R_base*2.0, scaleL*R_base*2.0, 1.0, 0.0, &axis1);
    d2.AddPlane (-17, Vec3_t(Center(0)-R_base*scaleL,Center(1),Center(2)), R, 2.0*L_z+height, scaleL*R_base*2.0, 1.0, ang, &axis1);
    d2.AddPlane (-18, Vec3_t(Center(0)+R_base*scaleL,Center(1),Center(2)), R, 2.0*L_z+height, scaleL*R_base*2.0, 1.0, ang, &axis1);
    d2.AddPlane (-19, Vec3_t(Center(0),Center(1)-R_base*scaleL,Center(2)), R, scaleL*R_base*2.0, 2.0*L_z+height, 1.0, M_PI/2.0, &axis0);
    d2.AddPlane (-20, Vec3_t(Center(0),Center(1)+R_base*scaleL,Center(2)), R, scaleL*R_base*2.0, 2.0*L_z+height, 1.0, M_PI/2.0, &axis0);
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
            d2.Particles[np]->Props.Gn = Gnb; // restitution coefficient
            d2.Particles[np]->Props.Mu = Muw; // frictional coefficient
        }
        if (d2.Particles[np]->Tag < -15)
        {
            d2.Particles[np]->v==Vec3_t(0.0,0.0,0.0);
      
            d2.Particles[np]->Props.Kn = Kn; // normal stiffness
            d2.Particles[np]->Props.Kt = Kt; // trangential stiffness
            d2.Particles[np]->Props.Gn = Gnw; // restitution coefficient
            d2.Particles[np]->Props.Mu = Muw; // frictional coefficient
        }
    }

    dtdem = ratiots*d2.CriticalDt(); //Calculating time step
    cout <<"Gnb"<< Gnb<<endl;
    cout <<"dtdem"<< dtdem<<endl;
    cout <<"ratiots"<< ratiots<<endl;
    cout <<"CACU"<< d2.CriticalDt()<<endl;
    cout <<"Nproc"<< Nproc<<endl;
    cout <<"Omp"<< omp_get_max_threads()<<endl;

    d2.Alpha = R/2.0; //Verlet distance
    d2.Solve(/*tf*/Tf2, dtdem, /*dtOut*/dtOut2, Setup, Report, "sandpile", 2, Nproc);
    d2.Save("Stage_cpu");

}
MECHSYS_CATCH
