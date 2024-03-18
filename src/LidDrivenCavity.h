#pragma once

#include <string>
using namespace std;

class SolverCG;


/**
 * @brief LidDrivenCavity is a class designed to simulate fluid flow within a two-dimensional
 * lid-driven cavity. This fluid problem is commonly used as a benchmark problem to evaluate CFD
 * algorithm performances.
*/
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
    void SetThreads(int Nt);

    void SetNeighbour(int rank_up, int rank_down, int rank_left, int rank_right);
    void DomainDecomposition();


    void Initialise();
    void Integrate();
    void WriteSolution(std::string file);
    void PrintConfiguration();

    // MPI Parallel processing
    void InitialiseParallel();
    


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
    void Advance(bool verbose_advance);

    // Parallel processing functions
    void IntegrateParallel();



private:

    double* v   = nullptr;
    double* s   = nullptr;

    double* u0   = nullptr;
    double* u1   = nullptr;

    double dt   = 0.01;
    double T    = 1.0;
    double dx;
    double dy;
    int    Nx   = 9;
    int    Ny   = 9;
    int    Npts = 81;
    int    Nx_inner;
    int    Ny_inner;
    int    Npts_inner;


    double Lx   = 1.0;
    double Ly   = 1.0;
    double Re   = 10;
    double U    = 1.0;
    double nu   = 0.1;

    // MPI Parallel processing
    bool parallel; // Serial or parallel processing
    int    Npts_local;
    int    Npts_local_buffer;
    int    Nx_remainder; // remainder of nodes in x direction
    int    Ny_remainder; // remainder of nodes in y direction
    int    Nx_local;
    int    Ny_local;

    int    offset_x;
    int    offset_y;

    double* v_local = nullptr; // local vorticity matrix
    double* s_local = nullptr; // local stream function matrix

    double* A_global_temp = nullptr; // temporary pointer storing global matrix


    bool verbose = true; // Display convergence and timestep detail during integration

    // MPI process rank and number of processes
    int root    = 0;
    int rank;
    int Nprocs;
    int Nprocs_sqrt; // Square root of nprocs for domain decomposition
    int rank_up, rank_down, rank_left, rank_right; // neighbour ranking


    double* buffer_up_send = nullptr;
    double* buffer_down_send = nullptr;
    double* buffer_up_recv = nullptr;
    double* buffer_down_recv = nullptr;

    double* buffer_left_send = nullptr;
    double* buffer_right_send = nullptr;
    double* buffer_left_recv = nullptr;
    double* buffer_right_recv = nullptr;


    SolverCG* cg = nullptr;

    void CleanUp();
    void UpdateDxDy();

    // Parallel processing functions
    void GatherDomain(double* A_local, double* A_global);
    void ScatterDomain(double* A_local, double* A_global);
    void DomainInterComunnication(double* A_local);


    void AdvanceParallel(bool verbose_advance);

    void ComputeBoundaryVorticityParallel();
    void ComputeInteriorVorticityParallel();
    void ComputeNextVorticityParallel();
    void ComputeLaplaceOperatorParallel();

    // Threading
    int Nt;


    void CreateU();
    void Printmatrix(int nx, int ny, double* A);
};

