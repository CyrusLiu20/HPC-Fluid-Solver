#pragma once

#include <string>
using namespace std;

class SolverCG;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetVerbose(bool verbose);


    void Initialise();
    void Integrate();
    void WriteSolution(std::string file);
    void PrintConfiguration();

    // Test cases (Developer use)
    double get_dx();
    double get_dy();
    double get_nu();
    int get_Npts();
    double* get_v();
    double* get_s();
    double* get_u0();
    double* get_u1();

    void IntegrateControl(double percentage);
    void Advance();




private:
    double* v   = nullptr;
    double* s   = nullptr;
    double* tmp = nullptr;

    double dt   = 0.01;
    double T    = 1.0;
    double dx;
    double dy;
    int    Nx   = 9;
    int    Ny   = 9;
    int    Npts = 81;
    double Lx   = 1.0;
    double Ly   = 1.0;
    double Re   = 10;
    double U    = 1.0;
    double nu   = 0.1;

    bool verbose = true;

    SolverCG* cg = nullptr;

    void CleanUp();
    void UpdateDxDy();
};

