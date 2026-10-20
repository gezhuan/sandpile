/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 ************************************************************************/

#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/linalg/matvec.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

using std::cout;
using std::endl;

struct HighOutputConfig
{
    bool enabled;
    double t1;
    double t2;
    double dt;

    HighOutputConfig() : enabled(false), t1(0.0), t2(0.0), dt(0.0) {}
};

struct InputData
{
    String ptype;
    bool Render;
    double L_x, L_y, L_z;
    double dt, Dp, R;
    double Tf1, dtOut1, Tf2, dtOut2;
    double height, height2;
    double R_base, rho_s, w_th, Muw;
    double Kn, Kt, Gn, Gnb;
    double Rp_base, Gnw, sizedis, MUP;
    double Ang, fraction, R_source;
    int Nump;
    double tp;
    int Anumper;
    String contactlaw;
    bool stage_packing, stage_avalanche;
    double scaleR;

    HighOutputConfig hfPacking;
    HighOutputConfig hfAvalanche;
};

struct UserData
{
    double height;
    double R;
    double rho_s;
    double acc;
    double tp;
    int ADDS;
    double Kn;
    double Kt;
    double Gn;
    double MUP;
    double timeperadd;
    double R_source;
    double w_th;
    double scaleR;
    String ptype;
    int Nump;
    int Anumper;
    size_t injectEverySteps;

    UserData()
        : height(0.0), R(0.0), rho_s(0.0), acc(0.0), tp(0.0), ADDS(0),
          Kn(0.0), Kt(0.0), Gn(0.0), MUP(0.0), timeperadd(0.0),
          R_source(0.0), w_th(0.0), scaleR(1.0), Nump(0), Anumper(1),
          injectEverySteps(1) {}
};

struct Circle
{
    double x, y, r;
};

static void ReadValue(std::ifstream &in, String &v) { in >> v; in.ignore(200,'\n'); }
static void ReadValue(std::ifstream &in, bool   &v) { in >> v; in.ignore(200,'\n'); }
static void ReadValue(std::ifstream &in, int    &v) { in >> v; in.ignore(200,'\n'); }
static void ReadValue(std::ifstream &in, double &v) { in >> v; in.ignore(200,'\n'); }

static InputData ReadInput(const String &filename)
{
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());

    std::ifstream in(filename.CStr());
    InputData p;

    ReadValue(in,p.ptype);
    ReadValue(in,p.Render);
    ReadValue(in,p.L_x);
    ReadValue(in,p.L_y);
    ReadValue(in,p.L_z);
    ReadValue(in,p.dt);
    ReadValue(in,p.Dp);
    ReadValue(in,p.R);
    ReadValue(in,p.Tf1);
    ReadValue(in,p.dtOut1);
    ReadValue(in,p.Tf2);
    ReadValue(in,p.dtOut2);
    ReadValue(in,p.height);
    ReadValue(in,p.height2);
    ReadValue(in,p.R_base);
    ReadValue(in,p.rho_s);
    ReadValue(in,p.w_th);
    ReadValue(in,p.Muw);
    ReadValue(in,p.Kn);
    ReadValue(in,p.Kt);
    ReadValue(in,p.Gn);
    ReadValue(in,p.Gnb);
    ReadValue(in,p.Rp_base);
    ReadValue(in,p.Gnw);
    ReadValue(in,p.sizedis);
    ReadValue(in,p.MUP);
    ReadValue(in,p.Ang);
    ReadValue(in,p.fraction);
    ReadValue(in,p.R_source);
    ReadValue(in,p.Nump);
    ReadValue(in,p.tp);
    ReadValue(in,p.Anumper);
    ReadValue(in,p.contactlaw);
    ReadValue(in,p.stage_packing);
    ReadValue(in,p.stage_avalanche);
    ReadValue(in,p.scaleR);
    ReadValue(in,p.hfPacking.enabled);
    ReadValue(in,p.hfPacking.t1);
    ReadValue(in,p.hfPacking.t2);
    ReadValue(in,p.hfPacking.dt);
    ReadValue(in,p.hfAvalanche.enabled);
    ReadValue(in,p.hfAvalanche.t1);
    ReadValue(in,p.hfAvalanche.t2);
    ReadValue(in,p.hfAvalanche.dt);

    return p;
}

static bool CirclesOverlap(const Circle &a, const Circle &b)
{
    const double dx = a.x-b.x;
    const double dy = a.y-b.y;
    const double rr = a.r+b.r;
    return dx*dx + dy*dy < rr*rr;
}

static std::vector<Circle> GenerateCircles(double bigR, double smallR, int N,
                                            int maxAttemptsPerCircle=1000,
                                            int maxTotalAttempts=10000)
{
    const unsigned seed = static_cast<unsigned>(
        std::chrono::system_clock::now().time_since_epoch().count());
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> angleDist(0.0,2.0*M_PI);
    std::uniform_real_distribution<double> radiusDist(0.0,1.0);

    std::vector<Circle> circles;
    circles.reserve(N);
    int totalAttempts = 0;

    for (int i=0; i<N && totalAttempts<maxTotalAttempts; ++i)
    {
        bool placed = false;
        for (int attempt=0; attempt<maxAttemptsPerCircle && totalAttempts<maxTotalAttempts;
             ++attempt, ++totalAttempts)
        {
            const double angle  = angleDist(generator);
            const double radius = std::sqrt(radiusDist(generator))*(bigR-smallR);
            Circle candidate = {radius*std::cos(angle), radius*std::sin(angle), smallR};

            bool valid = true;
            for (size_t j=0; j<circles.size(); ++j)
            {
                if (CirclesOverlap(candidate,circles[j]))
                {
                    valid = false;
                    break;
                }
            }

            if (valid)
            {
                circles.push_back(candidate);
                placed = true;
                break;
            }
        }

        if (!placed)
        {
            std::cerr << "Warning: only placed " << circles.size()
                      << " of " << N << " requested particles\n";
            break;
        }
    }

    return circles;
}

static void InitUserData(UserData &dat, const InputData &p, double height)
{
    dat.height   = height;
    dat.R        = p.R;
    dat.rho_s    = p.rho_s;
    dat.acc      = p.Dp;
    dat.tp       = p.tp;
    dat.Kn       = p.Kn;
    dat.Kt       = p.Kt;
    dat.Gn       = p.Gn;
    dat.MUP      = p.MUP;
    dat.R_source = p.R_source;
    dat.w_th     = p.w_th;
    dat.scaleR   = p.scaleR;
    dat.ptype    = p.ptype;
    dat.Nump     = p.Nump;
    dat.Anumper  = p.Anumper;

    if (p.dt<=0.0) throw new Fatal("dt must be positive");
    if (std::abs(dat.acc)<=0.0) throw new Fatal("Dp/acc must be non-zero");

    dat.timeperadd = std::sqrt(4.0*dat.R/std::abs(dat.acc));
    dat.injectEverySteps = static_cast<size_t>(
        std::ceil(dat.timeperadd*dat.tp/p.dt));
    if (dat.injectEverySteps<1) dat.injectEverySteps=1;
}

static void SetParticleMaterial(DEM::Particle *p, double Kn, double Kt,
                                double Gn, double Mu)
{
    p->Props.Kn = Kn;
    p->Props.Kt = Kt;
    p->Props.Gn = Gn;
    p->Props.Mu = Mu;
}

static void ConfigureParticles(DEM::Domain &d, const InputData &p)
{
    for (size_t i=0; i<d.Particles.Size(); ++i)
    {
        DEM::Particle *par = d.Particles[i];

        if (par->Tag==-1)
        {
            par->Ff = par->Props.m*Vec3_t(0.0,0.0,p.Dp);
            SetParticleMaterial(par,p.Kn,p.Kt,p.Gn,p.MUP);
        }
        else if (par->Tag==-15)
        {
            par->v = Vec3_t(0.0,0.0,0.0);
            par->FixVeloc();
            SetParticleMaterial(par,p.Kn,p.Kt,p.Gnb,p.Muw);
        }
        else if (par->Tag<-15)
        {
            par->v = Vec3_t(0.0,0.0,0.0);
            SetParticleMaterial(par,p.Kn,p.Kt,p.Gnw,p.Muw);
        }
    }
}

static void FixContainerPlanes(DEM::Domain &d)
{
    for (int tag=-16; tag>=-20; --tag) d.GetParticle(tag)->FixVeloc();
}

static void ConfigureHighOutput(DEM::Domain &d, const HighOutputConfig &hf,
                                const char *stageName)
{
    if (!hf.enabled) return;
    if (hf.t1<0.0 || hf.t2<hf.t1 || hf.dt<=0.0)
        throw new Fatal("Invalid %s high-frequency output window",stageName);

    // High-Frequency Standard Output V0.1:
    // reuse the DEM standard particle/contact writers inside Solve.
    d.SetHighFrequencyOutput(hf.t1,hf.t2,hf.dt);
    d.SetHighStandardOutput(true,true); // particles, contacts

    cout << stageName << " high-frequency standard output: ["
         << hf.t1 << ", " << hf.t2 << "] every " << hf.dt
         << " (particles + contacts)" << endl;
}

static void AddInjectedParticle(DEM::Domain &d, UserData &dat,
                                const Circle &c, Vec3_t &axis0)
{
    const Vec3_t pos(c.x,c.y,dat.height);

    if      (dat.ptype=="sphere") d.AddSphere(-1,pos,dat.R,dat.rho_s);
    else if (dat.ptype=="cube")   d.AddCube  (-1,pos,dat.w_th,dat.scaleR*dat.R,dat.rho_s);
    else if (dat.ptype=="tetra")  d.AddTetra (-1,pos,dat.w_th,dat.scaleR*dat.R,dat.rho_s);
    else if (dat.ptype=="rice")   d.AddRice  (-1,pos,dat.w_th,dat.scaleR*dat.R,dat.rho_s,M_PI/2.0,&axis0);
    else throw new Fatal("Unknown particle type <%s>",dat.ptype.CStr());

    const size_t idx = d.Particles.Size()-1;
    DEM::Particle *p = d.Particles[idx];
    SetParticleMaterial(p,dat.Kn,dat.Kt,dat.Gn,dat.MUP);
    p->Ff = p->Props.m*Vec3_t(0.0,0.0,dat.acc);
    p->v  = Vec3_t(0.0,0.0,0.0);
    p->InitializeVelocity(d.Dt);
    d.FreePar.Push(idx);
}

static void SetupInjection(DEM::Domain &d, void *UD)
{
    UserData &dat = *static_cast<UserData *>(UD);
    if (dat.ADDS>=dat.Nump || d.iter%dat.injectEverySteps!=0) return;

    const int requestN = std::min(dat.Anumper,dat.Nump-dat.ADDS);
    if (requestN<=0) return;

    std::vector<Circle> circles = GenerateCircles(dat.R_source*dat.R,dat.R,requestN);
    Vec3_t axis0(OrthoSys::e0);

    for (size_t i=0; i<circles.size(); ++i)
    {
        AddInjectedParticle(d,dat,circles[i],axis0);
        ++dat.ADDS;
    }

    if (!circles.empty()) d.UpdateContacts();
}

static size_t ContactLaw(const String &name)
{
    return name=="hertz" ? 1 : 0;
}

static void AddInitialParticle(DEM::Domain &d, const InputData &p,
                               double height, Vec3_t &axis0)
{
    const Vec3_t pos(0.0,0.0,height-16.0*p.R);
    if      (p.ptype=="sphere") d.AddSphere(-1,pos,p.R,p.rho_s);
    else if (p.ptype=="cube")   d.AddCube  (-1,pos,p.w_th,p.scaleR*p.R,p.rho_s);
    else if (p.ptype=="tetra")  d.AddTetra (-1,pos,p.w_th,p.scaleR*p.R,p.rho_s);
    else if (p.ptype=="rice")   d.AddRice  (-1,pos,p.w_th,p.scaleR*p.R,p.rho_s,M_PI/2.0,&axis0);
    else throw new Fatal("Unknown particle type <%s>",p.ptype.CStr());
}

static void BuildPackingGeometry(DEM::Domain &d, const InputData &p)
{
    const Vec3_t center(0.0,0.0,0.0);
    Vec3_t axis0(OrthoSys::e0);
    Vec3_t axis1(OrthoSys::e1);

    AddInitialParticle(d,p,p.height,axis0);

    const Vec3_t xmin(center(0)-1.2*p.R_base,center(1)-1.2*p.R_base,center(2));
    const Vec3_t xmax(center(0)+1.2*p.R_base,center(1)+1.2*p.R_base,center(2)+2.0*p.Rp_base);
    d.GenSpheresBox(-15,xmin,xmax,p.Rp_base,p.rho_s,"HCP",1234,1.0,p.sizedis);

    const double scaleL = 2.0;
    const double angle  = p.Ang*M_PI/180.0;
    d.AddPlane(-16,Vec3_t(0.0,0.0,-p.L_z),p.w_th,scaleL*p.R_base*2.0,scaleL*p.R_base*2.0,1.0,0.0,&axis1);
    d.AddPlane(-17,Vec3_t(-p.R_base*scaleL,0.0,0.0),p.w_th,2.0*p.L_z+p.height,scaleL*p.R_base*2.0,1.0,angle,&axis1);
    d.AddPlane(-18,Vec3_t( p.R_base*scaleL,0.0,0.0),p.w_th,2.0*p.L_z+p.height,scaleL*p.R_base*2.0,1.0,angle,&axis1);
    d.AddPlane(-19,Vec3_t(0.0,-p.R_base*scaleL,0.0),p.w_th,scaleL*p.R_base*2.0,2.0*p.L_z+p.height,1.0,M_PI/2.0,&axis0);
    d.AddPlane(-20,Vec3_t(0.0, p.R_base*scaleL,0.0),p.w_th,scaleL*p.R_base*2.0,2.0*p.L_z+p.height,1.0,M_PI/2.0,&axis0);
    FixContainerPlanes(d);

    for (size_t i=0; i<d.Particles.Size(); ++i)
    {
        DEM::Particle *par = d.Particles[i];
        if (par->Tag!=-15) continue;

        const double dx=par->x(0), dy=par->x(1), dz=par->x(2);
        if (dx*dx+dy*dy+dz*dz >= p.R_base*p.R_base) par->Tag=10;
    }

    Array<int> tags;
    tags.Push(10);
    d.DelParticles(tags);
}

static void RemoveParticlesBelowPackingBase(DEM::Domain &d, double R)
{
    for (size_t i=0; i<d.Particles.Size(); ++i)
    {
        if (d.Particles[i]->Tag==-1 && d.Particles[i]->x(2)<=-8.0*R)
            d.Particles[i]->Tag=11;
    }

    Array<int> tags;
    tags.Push(11);
    d.DelParticles(tags);
}

static void RunPacking(const InputData &p, size_t Nproc, size_t cl)
{
    UserData dat;
    InitUserData(dat,p,p.height);

    DEM::Domain d(&dat,cl);
    d.UserData=&dat;

    BuildPackingGeometry(d,p);
    ConfigureParticles(d,p);
    d.Alpha=p.R;
    ConfigureHighOutput(d,p.hfPacking,"Packing");

    cout << "Number of Particles " << d.Particles.Size() << endl;
    d.Solve(p.Tf1,p.dt,p.dtOut1,SetupInjection,NULL,"packing",2,Nproc);

    RemoveParticlesBelowPackingBase(d,p.R);
    d.Save("packing");
}

static void RunAvalanche(const InputData &p, size_t Nproc, size_t cl)
{
    UserData dat;
    InitUserData(dat,p,p.height2);

    DEM::Domain d(&dat,cl);
    d.UserData=&dat;
    d.Load("packing");

    FixContainerPlanes(d);
    ConfigureParticles(d,p);

    cout << "CACU "  << d.CriticalDt()      << endl;
    cout << "Nproc " << Nproc               << endl;
    cout << "Omp "   << omp_get_max_threads()<< endl;

    d.Alpha=p.R/2.0;
    ConfigureHighOutput(d,p.hfAvalanche,"Avalanche");
    d.Solve(p.Tf2,p.dt,p.dtOut2,SetupInjection,NULL,"sandpile",2,Nproc);
    d.Save("Stage_final");
}

int main(int argc, char **argv) try
{
    if (argc<2)
        throw new Fatal("This program must be called with one argument: the input file key without '.inp'.\nExample:\t%s filekey\n",argv[0]);

    const size_t Nproc = argc>=3 ? static_cast<size_t>(atoi(argv[2])) : 1;
    const String filekey(argv[1]);
    const InputData p = ReadInput(filekey+".inp");
    const size_t cl = ContactLaw(p.contactlaw);

    cout << "tp "              << p.tp              << endl;
    cout << "Rp_base "         << p.Rp_base         << endl;
    cout << "L_z "             << p.L_z             << endl;
    cout << "stage_packing "   << p.stage_packing   << endl;
    cout << "stage_avalanche " << p.stage_avalanche << endl;
    cout << "contactlaw "      << p.contactlaw      << endl;

    if (p.stage_packing)   RunPacking(p,Nproc,cl);
    if (p.stage_avalanche) RunAvalanche(p,Nproc,cl);
}
MECHSYS_CATCH