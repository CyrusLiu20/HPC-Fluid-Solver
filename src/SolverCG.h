#pragma once

/**
 * @brief SolverCG class for solving a linear system using the Conjugate Gradient method. It is used to
 * solve the laplace operator in the vorticity stream function eqution.
 * This class supports the usage for both serial and parallel computation.
 */
class SolverCG
{
public:
    SolverCG(int pNx, int pNy, double pdx, double pdy);
    ~SolverCG();

    // Conjugate gradient solver configuration
    void SetNeighbour(int rank_up, int rank_down, int rank_left, int rank_right);
    void SetOffset(int offset_x, int offset_y, int Nx_global, int Ny_global);

    // Find solution vector
    void Solve(double* b, double* x, bool verbose); // Legacy baseline code
    void SolveParallel(double* b, double* x, bool verbose);

private:

    // CG solver configuration parameters
    double dx;
    double dy;
    int Nx;
    int Ny;
    double* r;
    double* p;
    double* z;
    double* t;
    int iter_max; // Maximum iteration
    bool debug;

    // Parallel solving parameters
    int rank_up, rank_down, rank_left, rank_right; // Neighbour ranking
    int rank;
    int root = 0;
    int Nprocs;
    int i_start;
    int i_end;
    int j_start;
    int j_end;
    int offset_x;
    int offset_y;
    int Nx_global;
    int Ny_global;
    int Npts;
    double dx2i;
    double dy2i;
    double factor2;

    // Buffer memory for boundary values communication
    double* buffer_up_send = nullptr;
    double* buffer_down_send = nullptr;
    double* buffer_up_recv = nullptr;
    double* buffer_down_recv = nullptr;
    double* buffer_left_send = nullptr;
    double* buffer_right_send = nullptr;
    double* buffer_left_recv = nullptr;
    double* buffer_right_recv = nullptr;

    // Utilities
    void Printmatrix(int nx, int ny, double* A);
    void GatherDomain(double* A_local, double* A_global);

    // Serial computing
    void Precondition(double* p, double* t);
    void ImposeBC(double* p);

    // Parallel computing
    void DomainInterComunnication(double* A_local);
    void ApplyOperator(double* p, double* t);
    void PreconditionParallel(double* p, double* t);
    void ImposeBCParallel(double* p);
    double ComputeErrorGlobalParallel(unsigned int n, double *a);
    double ComputeDotGlobalParallel(unsigned int n, double *a, double *b);


};

