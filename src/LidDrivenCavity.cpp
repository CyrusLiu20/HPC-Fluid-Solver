#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <mpi.h>
using namespace std;

#include <cblas.h>

#define IDX(I,J) ((J)*Nx + (I))

#include "LidDrivenCavity.h"
#include "SolverCG.h"

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
    CleanUp();
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    this->Lx = xlen;
    this->Ly = ylen;
    UpdateDxDy();
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    this->Nx = nx;
    this->Ny = ny;

    // Interior matrix size
    this->Nx_inner = Nx - 2;
    this->Ny_inner = Ny - 2;
    this->Npts_inner = Nx_inner*Ny_inner;
    UpdateDxDy();
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    this->dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
    this->T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    this->Re = re;
    this->nu = 1.0/re;
}

void LidDrivenCavity::SetVerbose(bool verbose)
{
    this->verbose = verbose;
}

void LidDrivenCavity::InitialiseParallel()
{
    v_local = new double[Npts_local];
    s_local = new double[Npts_local];
    s = new double[Npts];
    // cg = new SolverCG(Nx_local,Ny_local,dx_local,dy_local);

    if(rank=root){
        v = new double[Npts];
    }

}

void LidDrivenCavity::DomainDecomposition()
{

    // Check if MPI is initialised
    int init_mpi;
    MPI_Initialized(&init_mpi);
    if (!init_mpi) {
        std::cerr << "Error: MPI not initialized" << std::endl;
        throw std::runtime_error("MPI not initialized");
    }
    else{
        // Retrieve current rank and process
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);

        parallel = Nprocs > 1;
        if(parallel){
            Nprocs_sqrt = sqrt(Nprocs);
            Nx_remainder = Nx_inner % Nprocs_sqrt;
            Ny_remainder = Ny_inner % Nprocs_sqrt;
            // std::cout << "Nx_inner: " << Nx_inner << std::endl;

            if(rank % Nprocs_sqrt < Nx_remainder){
                Nx_local = Nx_inner / Nprocs_sqrt + 1;
                offset_x = Nx_local * (rank % Nprocs_sqrt);
            }
            else{
                Nx_local = Nx_inner / Nprocs_sqrt;
                offset_x = Nx_local * (rank % Nprocs_sqrt) + Nx_remainder;
            }

            if(rank / Nprocs_sqrt < Ny_remainder){
                Ny_local = Ny_inner / Nprocs_sqrt + 1;
                offset_y = Ny_local * (rank / Nprocs_sqrt);
            }
            else{
                Ny_local = Ny_inner / Nprocs_sqrt;
                offset_y = Ny_local * (rank / Nprocs_sqrt) + Ny_remainder;
            }

            Npts_local = Nx_local*Ny_local;
            std::cout << "Rank: " << rank  << " | Nx_local: " << Nx_local << " | offset_x: " << offset_x << " | Ny_local: " << Ny_local << " | offset_y: " << offset_y << std::endl;

        }
    }
}

void LidDrivenCavity::ComputeBoundaryVorticityParallel(){

    if(rank=root){
        double dxi  = 1.0/dx;
        double dyi  = 1.0/dy;
        double dx2i = 1.0/dx/dx;
        double dy2i = 1.0/dy/dy;

        // Boundary node vorticity
        for (int i = 1; i < Nx-1; ++i) {
            // top
            v[IDX(i,0)]    = 2.0 * dy2i * (s[IDX(i,0)]    - s[IDX(i,1)]);
            // bottom
            v[IDX(i,Ny-1)] = 2.0 * dy2i * (s[IDX(i,Ny-1)] - s[IDX(i,Ny-2)])
                        - 2.0 * dyi*U;
        }
        for (int j = 1; j < Ny-1; ++j) {
            // left
            v[IDX(0,j)]    = 2.0 * dx2i * (s[IDX(0,j)]    - s[IDX(1,j)]);
            // right
            v[IDX(Nx-1,j)] = 2.0 * dx2i * (s[IDX(Nx-1,j)] - s[IDX(Nx-2,j)]);
        }
    }

}

void LidDrivenCavity::ComputeInteriorVorticityParallel(){

    if(parallel){
        double dxi  = 1.0/dx;
        double dyi  = 1.0/dy;
        double dx2i = 1.0/dx/dx;
        double dy2i = 1.0/dy/dy;

        for (int i = 0; i < Nx_local; ++i) {
            for (int j = 0; j < Ny_local; ++j) {
                v_local[IDX(i,j)] = dx2i*(
                        2.0 * s[IDX(i,j)] - s[IDX(i+1,j)] - s[IDX(i-1,j)])
                            + 1.0/dy/dy*(
                        2.0 * s[IDX(i,j)] - s[IDX(i,j+1)] - s[IDX(i,j-1)]);
            }
        }



    }

}

void LidDrivenCavity::Initialise()
{
    CleanUp();

    v   = new double[Npts]();
    s   = new double[Npts]();
    tmp = new double[Npts]();
    cg  = new SolverCG(Nx, Ny, dx, dy);
}

void LidDrivenCavity::Integrate()
{
    int NSteps = ceil(T/dt);
    for (int t = 0; t < NSteps; ++t)
    {
        if(verbose){
            std::cout << "Step: " << setw(8) << t
                    << "  Time: " << setw(8) << t*dt
                    << std::endl;
        }
        Advance();
    }
}

void LidDrivenCavity::WriteSolution(std::string file)
{
    double* u0 = new double[Nx*Ny]();
    double* u1 = new double[Nx*Ny]();
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u0[IDX(i,j)] =  (s[IDX(i,j+1)] - s[IDX(i,j)]) / dy;
            u1[IDX(i,j)] = -(s[IDX(i+1,j)] - s[IDX(i,j)]) / dx;
        }
    }
    for (int i = 0; i < Nx; ++i) {
        u0[IDX(i,Ny-1)] = U;
    }

    std::ofstream f(file.c_str());
    std::cout << "Writing file " << file << std::endl;
    int k = 0;
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            k = IDX(i, j);
            f << i * dx << " " << j * dy << " " << v[k] <<  " " << s[k] 
              << " " << u0[k] << " " << u1[k] << std::endl;
        }
        f << std::endl;
    }
    f.close();

    delete[] u0;
    delete[] u1;
}


void LidDrivenCavity::PrintConfiguration()
{
    cout << "Grid size: " << Nx << " x " << Ny << endl;
    cout << "Spacing:   " << dx << " x " << dy << endl;
    cout << "Length:    " << Lx << " x " << Ly << endl;
    cout << "Grid pts:  " << Npts << endl;
    cout << "Timestep:  " << dt << endl;
    cout << "Steps:     " << ceil(T/dt) << endl;
    cout << "Reynolds number: " << Re << endl;
    cout << "Linear solver: preconditioned conjugate gradient" << endl;
    cout << endl;
    if (nu * dt / dx / dy > 0.25) {
        cout << "ERROR: Time-step restriction not satisfied!" << endl;
        cout << "Maximum time-step is " << 0.25 * dx * dy / nu << endl;
        exit(-1);
    }
}


void LidDrivenCavity::CleanUp()
{
    if (v) {
        delete[] v;
        delete[] s;
        delete[] tmp;
        delete cg;
    }
}


void LidDrivenCavity::UpdateDxDy()
{
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    Npts = Nx * Ny;
}


void LidDrivenCavity::Advance()
{
    double dxi  = 1.0/dx;
    double dyi  = 1.0/dy;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;

    // Boundary node vorticity
    for (int i = 1; i < Nx-1; ++i) {
        // top
        v[IDX(i,0)]    = 2.0 * dy2i * (s[IDX(i,0)]    - s[IDX(i,1)]);
        // bottom
        v[IDX(i,Ny-1)] = 2.0 * dy2i * (s[IDX(i,Ny-1)] - s[IDX(i,Ny-2)])
                       - 2.0 * dyi*U;
    }
    for (int j = 1; j < Ny-1; ++j) {
        // left
        v[IDX(0,j)]    = 2.0 * dx2i * (s[IDX(0,j)]    - s[IDX(1,j)]);
        // right
        v[IDX(Nx-1,j)] = 2.0 * dx2i * (s[IDX(Nx-1,j)] - s[IDX(Nx-2,j)]);
    }

    // Compute interior vorticity
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            v[IDX(i,j)] = dx2i*(
                    2.0 * s[IDX(i,j)] - s[IDX(i+1,j)] - s[IDX(i-1,j)])
                        + 1.0/dy/dy*(
                    2.0 * s[IDX(i,j)] - s[IDX(i,j+1)] - s[IDX(i,j-1)]);
        }
    }

    // Time advance vorticity
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            v[IDX(i,j)] = v[IDX(i,j)] + dt*(
                ( (s[IDX(i+1,j)] - s[IDX(i-1,j)]) * 0.5 * dxi
                 *(v[IDX(i,j+1)] - v[IDX(i,j-1)]) * 0.5 * dyi)
              - ( (s[IDX(i,j+1)] - s[IDX(i,j-1)]) * 0.5 * dyi
                 *(v[IDX(i+1,j)] - v[IDX(i-1,j)]) * 0.5 * dxi)
              + nu * (v[IDX(i+1,j)] - 2.0 * v[IDX(i,j)] + v[IDX(i-1,j)])*dx2i
              + nu * (v[IDX(i,j+1)] - 2.0 * v[IDX(i,j)] + v[IDX(i,j-1)])*dy2i);
        }
    }

    // Sinusoidal test case with analytical solution, which can be used to test
    // the Poisson solver
    /*
    const int k = 3;
    const int l = 3;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            v[IDX(i,j)] = -M_PI * M_PI * (k * k + l * l)
                                       * sin(M_PI * k * i * dx)
                                       * sin(M_PI * l * j * dy);
        }
    }
    */

    // Solve Poisson problem
    cg->Solve(v, s, verbose);
}

double LidDrivenCavity::get_dx()
{
    return dx;
}

double LidDrivenCavity::get_dy()
{
    return dy;
}

int LidDrivenCavity::get_Npts()
{
    return Npts;
}

double LidDrivenCavity::get_nu()
{
    return nu;
}

double* LidDrivenCavity::get_v()
{
    return v;
}

double* LidDrivenCavity::get_s()
{
    return s;
}

// Extract u velocity matrix by ds/dy = u
double* LidDrivenCavity::get_u0()
{
    double* u0 = new double[Nx*Ny]();
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u0[IDX(i,j)] =  (s[IDX(i,j+1)] - s[IDX(i,j)]) / dy;
        }
    }
    // Boundary condition where u0 at the top row equals U
    for (int i = 0; i < Nx; ++i) {
        u0[IDX(i,Ny-1)] = U;
    }

    return u0;
}

// Extract v velocity matrix by ds/dx = -v
double* LidDrivenCavity::get_u1()
{
    double* u1 = new double[Nx*Ny]();
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u1[IDX(i,j)] = -(s[IDX(i+1,j)] - s[IDX(i,j)]) / dx;
        }
    }

    return u1;
}

// Integrate up to a certain time
void LidDrivenCavity::IntegrateControl(double percentage)
{
    int NSteps = ceil(T/dt*percentage);
    // std::cout << floor(T/dt*percentage) << std::endl;
    for (int t = 0; t < NSteps; ++t)
    {
        if(verbose){
            std::cout << "Step: " << setw(8) << t << "  Time: " << setw(8) << t*dt << std::endl;
        }
        Advance();
    }
}