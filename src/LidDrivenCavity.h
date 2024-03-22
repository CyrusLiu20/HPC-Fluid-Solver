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

    // Configuring fluid solver
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetVerbose(bool verbose);
    void SetThreads(int Nt);
    void SetNeighbour(int rank_up, int rank_down, int rank_left, int rank_right);

    // Utilities
    void WriteSolution(std::string file);
    void PrintConfiguration();

    // Serial computing
    void Initialise();
    void Integrate();

    // MPI Parallel processing
    void InitialiseParallel();
    void IntegrateParallel();
    void DomainDecomposition();

    // Test cases (Developer use)
    double get_dx();
    double get_dy();
    double get_nu();
    int get_Npts();
    double* get_v();
    double* get_s();
    double* get_u0();
    double* get_u1();
    void IntegrateControl(double percentage); // Integrating up to certain point
    void Advance(bool verbose_advance);

private:

    // Fluid solver configuration parameters
    double dt   = 0.01;
    double T    = 1.0;
    double dx;
    double dy;
    double Lx   = 1.0;
    double Ly   = 1.0;
    double Re   = 10;
    double U    = 1.0;
    double nu   = 0.1;
    int    Nx   = 9;
    int    Ny   = 9;
    int    Npts = 81;
    bool verbose = true; // Display convergence and timestep detail during integration

    SolverCG* cg = nullptr; // Conjugate gradient solver

    // Serial Computing
    double* v   = nullptr;
    double* v_new   = nullptr;
    double* s   = nullptr;

    // Utilities (Writing velocity to file)
    double* u0   = nullptr;
    double* u1   = nullptr;

    // MPI Parallel processing
    int    Npts_local;
    int    Nx_remainder; // remainder of nodes in x direction
    int    Ny_remainder; // remainder of nodes in y direction
    int    Nx_local;
    int    Ny_local;
    int    offset_x;
    int    offset_y;
    double* v_local = nullptr; // local vorticity matrix
    double* v_next_local = nullptr; // local vorticity matrix for next time step
    double* s_local = nullptr; // local stream function matrix
    double* A_global_temp = nullptr; // temporary pointer storing global matrix

    // MPI process rank and number of processes
    int root    = 0;
    int rank;
    int Nprocs;
    int Nprocs_sqrt; // Square root of nprocs for domain decomposition
    int rank_up, rank_down, rank_left, rank_right; // neighbour ranking

    // Threading
    int Nt;

    // Buffer memory for sharing boundary value
    double* buffer_up_send = nullptr;
    double* buffer_down_send = nullptr;
    double* buffer_up_recv = nullptr;
    double* buffer_down_recv = nullptr;

    double* buffer_left_send = nullptr;
    double* buffer_right_send = nullptr;
    double* buffer_left_recv = nullptr;
    double* buffer_right_recv = nullptr;

    // Utilities
    void CleanUp();
    void UpdateDxDy();
    void CreateU();
    void Printmatrix(int nx, int ny, double* A);

    // Parallel processing functions
    void GatherDomain(double* A_local, double* A_global);
    void ScatterDomain(double* A_local, double* A_global);
    void DomainInterComunnication(double* A_local);
    void AdvanceParallel(bool verbose_advance);
    void ComputeBoundaryVorticityParallel();
    void ComputeInteriorVorticityParallel();
    void ComputeNextVorticityParallel();
    void ComputeLaplaceOperatorParallel();

};

