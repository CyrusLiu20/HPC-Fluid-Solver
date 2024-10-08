/**
 * Lid-Driven Cavity problem serves as an important test case for validating computational models, understanding fundamental 
 * fluid dynamics principles, and gaining insights into complex flow phenomena in confined geometries.
 * 
 * Flowchart/Timeline:
 * 
 *  1-8:    Configuring LidDrivenCavity Class
 *  1   SetDomainSize()
 *  2   SetGridSize()
 *  3   SetTimeStep()
 *  4   SetReynoldsNumber()
 *  5   SetVerbose()
 *  6   SetThreads()
 *  7   SetNeighbour()
 *  8   DomainDecomposition()
 * 
 *  9-17:   Begin simulation
 *  9   PrintConfiguration()
 *  10  InitialiseParallel()
 *      11  AdvanceParallel()
 *          12 ComputeBoundaryVorticityParallel()
 *          13 ComputeInteriorVorticityParalle()
 *          14 DomainInterCommunication()
 *          15 ComputeNextVorticityParallel()
 *          16 DomainInterCommunication()
 *          17 ComputeLaplaceOperatorParallel() : Compute Laplace operator with SolverCg
 *      18  GatherDomain() : Assemble global voriticity and stream matrix for writing
 * 
 *  19-20   Outputing simulation results
 *  19  WriteSollution()
 *  20  CleanUp()
 *  
*/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <mpi.h>
#include <omp.h>

using namespace std;

#include <cblas.h>

/**
 * @brief Natural ordering based on global x and y coordinates
 */
#define IDX(I,J) ((J)*Nx + (I))
/**
 * @brief Natural ordering based on local x and y coordinates
 */
#define IDX_local(I_local,J_local) ((J_local)*Nx_local + (I_local))

#include "LidDrivenCavity.h"
#include "SolverCG.h"


/**
 * @brief Prints the content of a 2D matrix to the standard output.
 * @param nx The number of columns in the matrix.
 * @param ny The number of rows in the matrix.
 * @param A Pointer to the 1D array representing the matrix.
*/
void LidDrivenCavity::Printmatrix(int nx, int ny, double* A) {
    cout.precision(4);
    for (int j = ny-1; j >= 0; --j) {
        for (int i = 0; i < nx; ++i) {
            cout << setw(10) << A[j*nx+i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

/**
 * @brief Gathers the local domain data from each process and sends it to the root process to assemble the global domain.
 * The function ensures proper alignment and communication between processes in a parallel environment.
 * 
 * @param A_local Pointer to the local domain array (vorticity or stream function) of size Nx_local * Ny_local.
 * @param A_global Pointer to the global domain array (vorticity or stream function) at the root process of size Nx * Ny.
 * 
 * @details This function is responsible for gathering the local domain data from each process and sending it to the root process
 * to assemble the global domain. It takes into account the local offsets and grid sizes to ensure proper alignment and communication
 * between processes in a parallel environment
 * 
 * @see ScatterDomain()
 */
void LidDrivenCavity::GatherDomain(double* A_local, double* A_global){

    // Parameters for the (true) local domain without buffwe
    int Nx_local_temp = Nx_local;
    int Ny_local_temp = Ny_local;
    int i_start = 0;
    int i_end = Nx_local;
    int j_start = 0;
    int j_end = Ny_local;

    // Check if required to exclude buffer zone
    if (rank_up != -2) {Ny_local_temp--;j_end--;}
    if (rank_down != -2) {Ny_local_temp--;j_start++;}
    if (rank_left != -2) {Nx_local_temp--;i_start++;}
    if (rank_right != -2) {Nx_local_temp--;i_end--;}

    int Npts_local_temp = Nx_local_temp*Ny_local_temp;
    double* A_local_temp = new double[Npts_local_temp]; // Temporarily store local domain without buffer
    for(int i=0;i<Npts_local_temp;i++){A_local_temp[i]=0;} // Initialise to prevent memory issue

    // Extracting true local domain (without buffer)
    int k = 0;
    for(int j=j_start;j<j_end;j++){
        for(int i=i_start;i<i_end;i++){
            A_local_temp[k] = A_local[IDX_local(i,j)];
            k++;
        }
    }

    if(rank!=root){
        // All non root processes send the number of parameters then the true local domain
        int send_buf[4] = {Nx_local_temp, Ny_local_temp, offset_x, offset_y};
        MPI_Send(send_buf,4,MPI_INT,root,0,MPI_COMM_WORLD);
        MPI_Send(A_local_temp, Npts_local_temp, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
    }
    else{
        int index_global, index_local;
        for (int src = 0; src < Nprocs; ++src) {

            // Unpacking number of elements to receive and the true local domain from non-root processes
            if (src != root) {
                int recv_buf[4];
                MPI_Recv(recv_buf, 4, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Number of elements to receive
                double* A_local_temp_recv = new double[recv_buf[0]*recv_buf[1]];

                // Receiving true local domain from all processes
                MPI_Recv(A_local_temp_recv, recv_buf[0]*recv_buf[1], MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int i=0;i<recv_buf[0];i++){
                    for(int j=0;j<recv_buf[1];j++){
                        index_global = (j + recv_buf[3])*Nx + (i + recv_buf[2]); // Global coordinates relative to local domain
                        index_local = (j)*recv_buf[0] + (i); // Local coordinates of the local domain
                        A_global[index_global] = A_local_temp_recv[index_local];
                    }
                }
                delete[] A_local_temp_recv;
            }
            else{
                // Assembling all true local domain based on local relative coordinates to A_global
                for(int i=i_start;i<i_end;i++){
                    for(int j=j_start;j<j_end;j++){
                        index_global = (j + offset_y)*Nx + (i + offset_x);
                        A_global[index_global] = A_local[IDX_local(i,j)];
                    
                    }
                }
            }
                
        }
    }
    // Synchronize all processes
    MPI_Barrier(MPI_COMM_WORLD);

    delete[] A_local_temp;
}

LidDrivenCavity::LidDrivenCavity()
{
}

/**
 * @brief Delete unnecessary memory for good practices
*/
LidDrivenCavity::~LidDrivenCavity()
{
    CleanUp();
}

/**
 * @brief Configure global domain size
 * @param xlen (double) global domain length in x-direction
 * @param ylen (double) global domain length in y-direction
*/
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    this->Lx = xlen;
    this->Ly = ylen;
    UpdateDxDy();
}

/**
 * @brief Configure global grid size and update mesh
 * @param nx (double) number of global grid points in x-direction
 * @param ny (double) number of global grid points in y-direction
*/
void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    this->Nx = nx;
    this->Ny = ny;
    UpdateDxDy();
}

/**
 * @brief Configure grid spacing based on domain size and number of grid points
*/
void LidDrivenCavity::UpdateDxDy()
{
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    Npts = Nx * Ny;
}


/**
 * @brief Configure time step
 * @param deltat (double) simulation time step
*/
void LidDrivenCavity::SetTimeStep(double deltat)
{
    this->dt = deltat;
}

/**
 * @brief Configure total simulation time
 * @param finalt (double) total simulation time
*/
void LidDrivenCavity::SetFinalTime(double finalt)
{
    this->T = finalt;
}

/**
 * @brief Configure simulation reynolds number
 * @param re (double) reynolds number
*/
void LidDrivenCavity::SetReynoldsNumber(double re)
{
    this->Re = re;
    this->nu = 1.0/re;
}

/**
 * @brief Configure verbosity option to display simulation time step and CG convergence detail
 * @param verbose (bool) verbosity option
*/
void LidDrivenCavity::SetVerbose(bool verbose)
{
    this->verbose = verbose;
}

/**
 * @brief Configure the number of threads
 * @param Nt (int) number of threads
*/
void LidDrivenCavity::SetThreads(int Nt)
{
    this->Nt = Nt;
}

/**
 * @brief Configures the neighboring process ranks for inter-process communication.
 * 
 * This function sets the ranks of the neighboring processes (up, down, left, and right) 
 * for the current process to facilitate inter-process communication in a distributed environment.
 * 
 * @param rank_up    Rank of the neighboring process above the current process.
 * @param rank_down  Rank of the neighboring process below the current process.
 * @param rank_left  Rank of the neighboring process to the left of the current process.
 * @param rank_right Rank of the neighboring process to the right of the current process.
 */
void LidDrivenCavity::SetNeighbour(int rank_up, int rank_down, int rank_left, int rank_right)
{
    this->rank_up = rank_up;
    this->rank_down = rank_down;
    this->rank_left = rank_left;
    this->rank_right = rank_right;

}

/**
 * @brief Performs domain decomposition for parallel computation and ensures buffer zone for information sharing
 * 
 * This function divides the computational domain into smaller subdomains for parallel 
 * computation using the MPI library. Each process is assigned a portion of the domain 
 * based on its rank, allowing for parallel execution of the Lid-Driven Cavity simulation.
 * 
 * The domain decomposition ensures that each process works on a distinct portion of 
 * the domain, facilitating efficient parallelization of the computation.
 */
void LidDrivenCavity::DomainDecomposition()
{

    // Retrieve current rank and process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);

    // Calculate remainder node points
    Nprocs_sqrt = sqrt(Nprocs);
    Nx_remainder = Nx % Nprocs_sqrt;
    Ny_remainder = Ny % Nprocs_sqrt;

    // Determine the local domain size and offset for the current process in x-direction
    if(rank % Nprocs_sqrt < Nx_remainder){
        Nx_local = Nx / Nprocs_sqrt + 1;
        offset_x = Nx_local * (rank % Nprocs_sqrt);
    }
    else{
        Nx_local = Nx / Nprocs_sqrt;
        offset_x = Nx_local * (rank % Nprocs_sqrt) + Nx_remainder;
    }

    // Determine the local domain size and offset for the current process in x-direction
    if(rank / Nprocs_sqrt < Ny_remainder){
        Ny_local = Ny / Nprocs_sqrt + 1;
        offset_y = Ny_local * (rank / Nprocs_sqrt);
    }
    else{
        Ny_local = Ny / Nprocs_sqrt;
        offset_y = Ny_local * (rank / Nprocs_sqrt) + Ny_remainder ;
    }

    // Buffer zone to account for shared rows/columns of neighbouring ranks
    if (rank_up != -2) {Ny_local++;}
    if (rank_down != -2) {Ny_local++;}
    if (rank_left != -2) {Nx_local++;}
    if (rank_right != -2) {Nx_local++;}

    Npts_local = Nx_local*Ny_local; // Number of elements in local domain
}

/**
 * @brief Initializes the LidDrivenCavity solver for parallel execution.
 * 
 * Allocates memory for local vorticity and stream function arrays, as well as buffers required for inter-domain communication.
 * Initializes all arrays to prevent arbitrary and unpredictable values.
 */
void LidDrivenCavity::InitialiseParallel()
{
    // Allocating memory for local vorticity and stream function
    v_local = new double[Npts_local];
    v_next_local = new double[Npts_local];
    s_local = new double[Npts_local];

    // Initialise conjugate gradient solver
    cg = new SolverCG(Nx_local,Ny_local,dx,dy);
    cg->SetNeighbour(rank_up,rank_down,rank_left,rank_right);
    cg->SetOffset(offset_x,offset_y,Nx,Ny);

    // Initialize elements for safety
    for (int i=0;i<Npts_local;i++) {
        v_local[i] = 0;
        s_local[i] = 0;
    }    

    // Initialise array to prevent memory issue
    if(rank==root){
        u0 = new double[Npts];
        u1 = new double[Npts];
        for (int i=0;i<Npts;i++) {u0[i] = 0;u1[i] = 0;}
    }
    if(rank==root){
        // Purely for gathering domain purposes and writing solution
        v = new double[Npts];
        s = new double[Npts];
        for (int i=0;i<Npts;i++) {v[i] = 0;s[i]=0;}
    }

    // Buffer for send/receive boundary information for domain inter-communication (Top and Bottom)
    buffer_up_send = new double[Nx_local];
    buffer_down_send = new double[Nx_local];
    buffer_up_recv = new double[Nx_local];
    buffer_down_recv = new double[Nx_local];
    for (int i=0;i<Nx_local;i++) {
        buffer_up_send[i] = 0;
        buffer_down_send[i] = 0;
        buffer_up_recv[i] = 0;
        buffer_down_recv[i] = 0;
    }    

    // Buffer for send/receive boundary information for domain inter-communication (Left and Right)
    buffer_left_send = new double[Ny_local];
    buffer_right_send = new double[Ny_local];
    buffer_left_recv = new double[Ny_local];
    buffer_right_recv = new double[Ny_local];
    for (int j=0;j<Ny_local;j++) {
        buffer_left_send[j] = 0;
        buffer_right_send[j] = 0;
        buffer_left_recv[j] = 0;
        buffer_right_recv[j] = 0;
    }    
}

/**
 * @brief Allocate memory for global u and v velocity and initialising to prevent memory leakage
*/
void LidDrivenCavity::CreateU(){
    u0 = new double[Npts];
    u1 = new double[Npts];
    for (int i=0;i<Npts;i++) {u0[i] = 0;u1[i] = 0;}

}

/**
 * @brief Integrates the Lid-Driven Cavity simulation over time in parallel.
 * This function iterates over time steps, advancing the simulation state at each step.
 * Gathers all local vorticity and stream for writing purposes
 * 
 * @details
 * The simulation proceeds for a total of `NSteps` time steps, where `NSteps` is determined
 * based on the total simulation time (`T`) and the time step size (`dt`).
 * 
 * Finally, the local simulation domain data is gathered from each process to reconstruct the 
 * complete simulation state for post-processing and analysis.
 * 
 *  @see AdvanceParallel(), ComputeBoundaryVorticityParallel(), ComputeInteriorVorticityParallel(), ComputeNextVorticityParallel(), ComputeLaplaceOperatorParallel()
 */
void LidDrivenCavity::IntegrateParallel(){
    int NSteps = ceil(T/dt);

    bool verbose_advance = false;
    for (int t = 0; t < NSteps; ++t)
    {
        if(verbose && rank==root){
            std::cout << "Step: " << setw(8) << t
                    << "  Time: " << setw(8) << t*dt
                    << std::endl;
        }

        // For debugging purposes
        if(t==-1){verbose_advance=true;}
        else{verbose_advance=false;}
        AdvanceParallel(verbose_advance); // Compute vorticity and steam function
    }

    // Assembling all local vorticity and stream function at root for writing purposes
    GatherDomain(v_local,v);
    GatherDomain(s_local,s);
}



/**
 * @brief Advances the Lid-Driven Cavity simulation to the next time step in parallel.
 * 
 * This function performs the necessary computations to advance the Lid-Driven Cavity 
 * simulation to the next time step in a parallel environment. It includes several steps 
 * such as computing boundary and interior vorticity, domain inter-communication, computing 
 * the Laplace operator, and updating the vorticity and stream function matrices.
 * 
 *  @see IntegrateParallel(), ComputeBoundaryVorticityParallel(), ComputeInteriorVorticityParallel(), ComputeNextVorticityParallel(), ComputeLaplaceOperatorParallel()
*/
void LidDrivenCavity::AdvanceParallel(bool verbose_advance){

    // Computing boundary and interior vorticity in parallel mode
    ComputeBoundaryVorticityParallel();
    ComputeInteriorVorticityParallel();

    if(rank==root&&verbose_advance){
        std::cout << "(Interior vorticity) vorticity matrix" << std::endl;
        Printmatrix(Nx, Ny, v);
    }
    DomainInterComunnication(v_local); // Sharing boundary values

    // Computing vorticity at next time step
    ComputeNextVorticityParallel();
    DomainInterComunnication(v_local); //  Sharing boundary values

    if(rank==root&&verbose_advance){
        std::cout << "(Next vorticity) vorticity matrix" << std::endl;
        Printmatrix(Nx, Ny, v);
        std::cout << "(Next vorticity) stream function matrix" << std::endl;
        Printmatrix(Nx, Ny, s);
    } 

    // Solution of the Poisson problem for stream function at next time step
    ComputeLaplaceOperatorParallel();

    if(rank==root&&verbose_advance){
        std::cout << "(Laplace operator) vorticity matrix" << std::endl;
        Printmatrix(Nx, Ny, v);
        std::cout << "(Laplace operator) stream function matrix" << std::endl;
        Printmatrix(Nx, Ny, s);
    }   

}


/**
 * @brief Computes the boundary vorticity for the Lid-Driven Cavity simulation in parallel.
 * 
 * This function calculates the boundary vorticity values of the local domain in a parallel 
 * environment. It computes the vorticity values along the top, bottom, left, and right 
 * boundaries of the local domain based on the given stream function values.
 * 
 * @details The boundary vorticity is computed using finite difference approximations 
 * and is updated directly in the local vorticity array. Boundary conditions are applied 
 * according to the coordinate of current process.
*/
void LidDrivenCavity::ComputeBoundaryVorticityParallel()
{

    double dyi  = 1.0/dy;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;

    // Top boundary values
    if (rank_up==-2){
        for (int i=1;i<Nx_local-1;++i){
            v_local[IDX_local(i,Ny_local-1)] = 2.0 * dy2i * (s_local[IDX_local(i,Ny_local-1)] - s_local[IDX_local(i,Ny_local-2)])- 2.0 * dyi*U;
        }
    }

    // Bottom boundary values
    if (rank_down==-2){
        for (int i=1;i<Nx_local-1;++i){
            v_local[IDX_local(i,0)] = 2.0 * dy2i * (s_local[IDX_local(i,0)] - s_local[IDX_local(i,1)]);
        }
    }

    // Left boundary values
    if (rank_left==-2){
        for (int j=1;j<Ny_local-1;++j){
            v_local[IDX_local(0,j)]    = 2.0 * dx2i * (s_local[IDX_local(0,j)] - s_local[IDX_local(1,j)]);
        }
    }

    // Right boundary values
    if (rank_right==-2){
        for (int j=1;j<Ny_local-1;++j){
            v_local[IDX_local(Nx_local-1,j)] = 2.0 * dx2i * (s_local[IDX_local(Nx_local-1,j)]    - s_local[IDX_local(Nx_local-2,j)]);
        }
    }
}


/**
 * @brief Computes the interior vorticity for the Lid-Driven Cavity simulation in parallel.
 * 
 * This function calculates the interior vorticity values of the local domain in a parallel 
 * environment. It computes the vorticity values for interior grid points based on the 
 * given stream function values using finite difference approximations.
 * 
 * @details The interior vorticity values are updated directly in the local vorticity array 
 * based on the computed finite difference equations. The computation is performed for 
 * interior grid points excluding boundary points to avoid redundant computation.
 * 
*/
void LidDrivenCavity::ComputeInteriorVorticityParallel(){
    double dx2i = 1.0/dx/dx;
    for (int i=1;i<Nx_local-1;i++){
        for (int j=1;j<Ny_local-1;j++){
            v_local[IDX_local(i,j)] = dx2i*(
                    2.0 * s_local[IDX_local(i,j)] - s_local[IDX_local(i+1,j)] - s_local[IDX_local(i-1,j)])
                        + 1.0/dy/dy*(
                    2.0 * s_local[IDX_local(i,j)] - s_local[IDX_local(i,j+1)] - s_local[IDX_local(i,j-1)]);
        }
    }
}



/**
 * @brief Computes the next vorticity values for the Lid-Driven Cavity simulation in parallel.
 * 
 * This function advances the vorticity values of the local domain to the next time step 
 * using the given stream function values and the current vorticity values. It implements 
 * the time integration scheme for the vorticity equation using finite difference 
 * approximations and updates the vorticity values accordingly.
*/
void LidDrivenCavity::ComputeNextVorticityParallel(){
    double dxi  = 1.0/dx;
    double dyi  = 1.0/dy;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    // Time advance vorticity
    for (int i=1;i<Nx_local-1;i++){
        for (int j=1;j<Ny_local-1;j++){ // Check whether to use another variable or not
            // v_local[IDX_local(i,j)] = v_local[IDX_local(i,j)] + dt*(
            v_next_local[IDX_local(i,j)] = v_local[IDX_local(i,j)] + dt*(
                ( (s_local[IDX_local(i+1,j)] - s_local[IDX_local(i-1,j)]) * 0.5 * dxi
                 *(v_local[IDX_local(i,j+1)] - v_local[IDX_local(i,j-1)]) * 0.5 * dyi)
              - ( (s_local[IDX_local(i,j+1)] - s_local[IDX_local(i,j-1)]) * 0.5 * dyi
                 *(v_local[IDX_local(i+1,j)] - v_local[IDX_local(i-1,j)]) * 0.5 * dxi)
              + nu * (v_local[IDX_local(i+1,j)] - 2.0 * v_local[IDX_local(i,j)] + v_local[IDX_local(i-1,j)])*dx2i
              + nu * (v_local[IDX_local(i,j+1)] - 2.0 * v_local[IDX_local(i,j)] + v_local[IDX_local(i,j-1)])*dy2i);
        }
    }
}


/**
 * @brief Computes the Laplace operator for the Lid-Driven Cavity simulation in parallel.
 * 
 * This function computes the Laplace operator (Laplacian) of the stream function 
 * in the local domain using a parallel Conjugate Gradient (CG) solver. It solves 
 * the Poisson equation ∇^2(s) = -v in the local domain, where ∇^2 represents the 
 * Laplacian operator, s is the stream function, and v is the vorticity
*/
void LidDrivenCavity::ComputeLaplaceOperatorParallel()
{
    // Solve Poisson problem
    cg->SolveParallel(v_next_local, s_local, verbose);

}


/**
 * @brief Manages inter-domain communication for the Lid-Driven Cavity simulation.
 * 
 * This function handles communication between neighboring processes in the parallel 
 * domain decomposition of the Lid-Driven Cavity simulation. It exchanges boundary 
 * information between adjacent domains to ensure consistency in the simulation results.
 * 
 * @param A_local Pointer to the local domain array containing the data to be communicated.
 * 
 * @see DomainDecomposition()
*/
void LidDrivenCavity::DomainInterComunnication(double* A_local){

    // Top and bottom communication
    // Sending upper edge of (true) local domain if neighbour exist
    if(rank_up!=-2){
        int j_second_upper = (Ny_local-1)-1;
        for(int i=0;i<Nx_local;i++){
            buffer_up_send[i] = A_local[IDX_local(i,j_second_upper)];
        }
    }
    // Sending bottom edge of (true) local domain if neighbour exist
    if(rank_down!=-2){
        int j_second_lower = (0)+1;
        for(int i=0;i<Nx_local;i++){
            buffer_down_send[i] = A_local[IDX_local(i,j_second_lower)];
        }
    }

    // Sending buffer data with neighboring processes if exist
    if(rank_up!=-2){
        MPI_Send(buffer_up_send, Nx_local, MPI_DOUBLE, rank_up, 0, MPI_COMM_WORLD);
    }
    if(rank_down!=-2){
        MPI_Send(buffer_down_send, Nx_local, MPI_DOUBLE, rank_down, 0, MPI_COMM_WORLD);
    }

    // Receiving buffer data with neighboring processes if exist
    if(rank_down!=-2){
        MPI_Recv(buffer_down_recv, Nx_local, MPI_DOUBLE, rank_down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(rank_up!=-2){
        MPI_Recv(buffer_up_recv, Nx_local, MPI_DOUBLE, rank_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Synchronise all processes
    MPI_Barrier(MPI_COMM_WORLD);

    // Unpacking all received data
    for(int i=0;i<Nx_local;i++){
        if (rank_up != -2) {A_local[IDX_local(i,Ny_local-1)] = buffer_up_recv[i];}
        if (rank_down != -2) {A_local[IDX_local(i,0)] = buffer_down_recv[i];}
    }

    // Left and right communication
    // Sending left edge of (true) local domain if neighbour exist
    if(rank_left!=-2){
        int i_second_left = (0) + 1;
        for(int j=0;j<Ny_local;j++){
            buffer_left_send[j] = A_local[IDX_local(i_second_left,j)];
        }
    }
    // Sending right edge of (true) local domain if neighbour exist
    if(rank_right!=-2){
        int i_second_right = (Nx_local-1) - 1;
        for(int j=0;j<Ny_local;j++){
            buffer_right_send[j] = A_local[IDX_local(i_second_right,j)];
        }
    }

    // Sending buffer data with neighboring processes if exist
    if(rank_left!=-2){
        MPI_Send(buffer_left_send, Ny_local, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD);
    }
    if(rank_right!=-2){
        MPI_Send(buffer_right_send, Ny_local, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD);
    }

    // Receiving buffer data with neighboring processes if exist
    if(rank_right!=-2){
        MPI_Recv(buffer_right_recv, Ny_local, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(rank_left!=-2){
        MPI_Recv(buffer_left_recv, Ny_local, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
    // Synchronise all processes
    MPI_Barrier(MPI_COMM_WORLD);

    // Unpacking all received data
    for(int j=0;j<Ny_local;j++){
        if (rank_left != -2) {A_local[IDX_local(0,j)] = buffer_left_recv[j];}
        if (rank_right != -2) {A_local[IDX_local(Nx_local-1,j)] = buffer_right_recv[j];}
    }       
}


/**
 * @brief Writes the solution data to a file.
 * 
 * This function writes the solution data, including vorticity, stream function, 
 * and velocity components, to a specified file in a structured format. It calculates 
 * the velocity components (u and v) from the stream function using finite differences 
 * and writes all solution data to the file.
 * 
 * @param file The name of the file to which the solution data will be written.
 * 
 * @note This function assumes that the solution data arrays (v, s, u0, u1) have been 
 * properly initialized and filled with appropriate values prior to the function call.
 * 
*/
void LidDrivenCavity::WriteSolution(std::string file)
{

    for(int i=0;i<Npts;i++){u0[i] = 0; u1[i] = 0;}
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u0[IDX(i,j)] =  (s[IDX(i,j+1)] - s[IDX(i,j)]) / dy; // Computing U-velocity from stream
            u1[IDX(i,j)] = -(s[IDX(i+1,j)] - s[IDX(i,j)]) / dx; // Computing V-velocity from stream
        }
    }

    for (int i = 0; i < Nx; ++i) {
        u0[IDX(i,Ny-1)] = U; // Boundary condition due to moving top lid
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
}


/**
 * @brief Prints the configuration parameters of the Lid-Driven Cavity simulation.
 * 
 * This function prints various configuration parameters of the Lid-Driven Cavity simulation,
 * including grid size, spacing, length, number of grid points, timestep, number of steps,
 * Reynolds number, number of processes, number of threads, and the linear solver used.
 * 
 * Additionally, it checks if the time-step restriction is satisfied based on the Courant-Friedrichs-Lewy (CFL) condition,
 * and displays an error message if the condition is violated, indicating that the maximum allowable time-step has been exceeded.
*/
void LidDrivenCavity::PrintConfiguration()
{
    cout << "Grid size: " << Nx << " x " << Ny << endl;
    cout << "Spacing:   " << dx << " x " << dy << endl;
    cout << "Length:    " << Lx << " x " << Ly << endl;
    cout << "Grid pts:  " << Npts << endl;
    cout << "Timestep:  " << dt << endl;
    cout << "Steps:     " << ceil(T/dt) << endl;
    cout << "Reynolds number: " << Re << endl;
    cout << "Number of processes: " << Nprocs << endl;
    cout << "Number of threads: " << Nt << endl;
    cout << "Linear solver: preconditioned conjugate gradient" << endl;

    cout << endl;
    if (nu * dt / dx / dy > 0.25) {
        cout << "ERROR: Time-step restriction not satisfied!" << endl;
        cout << "Maximum time-step is " << 0.25 * dx * dy / nu << endl;
        exit(-1);
    }
}


/**
 * @brief Cleans up memory allocated for vorticity, stream function, communication buffer, and the conjugate gradient solver.
*/
void LidDrivenCavity::CleanUp()
{
    if (v) {delete[] v;}
    if (s) {delete[] s;}
    if (v_new) {delete[] v_new;}

    if (v_local) {delete[] v_local;}
    if (s_local) {delete[] s_local;}
    if (v_next_local) {delete[] v_next_local;}
    if (cg) {delete cg;}

    if(buffer_up_send){
        delete[] buffer_up_send;
        delete[] buffer_down_send;
        delete[] buffer_up_recv;
        delete[] buffer_down_recv;

        delete[] buffer_left_send;
        delete[] buffer_right_send;
        delete[] buffer_left_recv;
        delete[] buffer_right_recv;
    }

    if(u0){delete[] u0;}
    if(u1){delete[] u1;}
}


/**
 * @brief (Primarily for debugging purposes, not used in the code) Scatters the global domain data from the root process to all other processes, ensuring proper distribution of data.
 * 
 * The function scatters the global matrix to all other processes.
 * It ensures proper alignment and distribution of data based on the process rank and offset values.
 * 
 * @param A_local Pointer to the local domain array (vorticity or stream function) of size Nx_local * Ny_local.
 * @param A_global Pointer to the global domain array (vorticity or stream function) at the root process of size Nx * Ny.
 * 
 * @deprecated This function is deprecated and will be removed in future versions.
 * 
 * @see GatherDomain()
 */
void LidDrivenCavity::ScatterDomain(double* A_local, double* A_global)
{

    if(rank==root){
        for (int src = 1; src < Nprocs; ++src) {
            MPI_Send(A_global, Npts, MPI_DOUBLE, src, 0, MPI_COMM_WORLD); // Sending A_global to each rank from root
        }
        int index_global, index_local;
        for (int i=0;i<Nx_local;i++){
            for (int j=0;j<Ny_local;j++){
                index_global = (j + 0)*Nx + (i + 0); // Global equivalent coordinate relative to local domain
                index_local = IDX_local(i,j); // Local coordinates   
                A_local[index_local] = A_global[index_global];
            }
        }
    } 
    else{
        double* A_global_temp = new double[Npts]; // Store global domain temporarily
        int index_global, index_local;
        MPI_Recv(A_global_temp, Npts, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Receive global domain from root

        // Relative position of local origin without buffer
        int offset_y_temp = offset_y;
        int offset_x_temp = offset_x;
        if(rank_down!=-2){offset_y_temp--;}
        if(rank_left!=-2){offset_x_temp--;}

        // Unpacking received global domain
        for (int i=0;i<Nx_local;i++){
            for (int j=0;j<Ny_local;j++){
                index_global = (j + offset_y_temp)*Nx + (i + offset_x_temp); // Global equivalent coordinate relative to local domain
                index_local = IDX_local(i,j); // Local coordinates           
                A_local[index_local] = A_global_temp[index_global];
            }
        }
        delete[] A_global_temp;
    }

}




















































//////////////////////////////////////////    Legacy baseline code    //////////////////////////////////////////


/**
 * @brief Initializes the LidDrivenCavity object by allocating memory for vorticity, stream function, and creating a solver object.
 * Also initializes velocity components.
 * 
 * @details This function cleans up any existing memory, allocates memory for vorticity (`v`), stream function (`s`), and creates a
 * conjugate gradient solver object (`cg`). Additionally, it initializes velocity components using the CreateU() function.
 * 
 * @deprecated This function is deprecated and will be removed in future versions. Consider using the InitialiseParallel() function for parallel initialization.
 * 
 * @see CleanUp(), CreateU(), InitialiseParallel()
 */
void LidDrivenCavity::Initialise()
{
    CleanUp();

    v   = new double[Npts]();
    v_new = new double[Npts]();
    s   = new double[Npts]();
    cg  = new SolverCG(Nx, Ny, dx, dy);
    CreateU();
}


/**
 * @brief Integrates the LidDrivenCavity object over time using a simple time-stepping scheme.
 * 
 * @details This function advances the simulation over time by iterating through a specified number of time steps.
 * Each time step involves calling the Advance function to update the solution. If the verbose mode is enabled,
 * it prints the step number and corresponding time at each iteration.
 * 
 * @deprecated This function is deprecated and will be removed in future versions. Consider using the IntegrateParallel() function for parallel integration.
 * 
 * @see Advance(), IntegrateParallel()
 */
void LidDrivenCavity::Integrate()
{
    int NSteps = ceil(T/dt);
    bool verbose_advance = false;
    for (int t = 0; t < NSteps; ++t)
    {
        if(verbose){
            std::cout << "Step: " << setw(8) << t
                    << "  Time: " << setw(8) << t*dt
                    << std::endl;
        }

        if(t==-1){verbose_advance=true;}
        else{verbose_advance=false;}
        Advance(verbose_advance);
    }
}

/**
 * @brief Integrates the Lid-Driven Cavity simulation over time in serial.
 * This function iterates over time steps, advancing the simulation state at each step.
 * 
 * @details
 * The simulation proceeds for a total of `NSteps` time steps, where `NSteps` is determined
 * based on the total simulation time (`T`) and the time step size (`dt`).
 * 
 * During each time step, the simulation state is advanced by calling the `Advance` function.
 * Optionally, verbose output can be enabled for the first time step for detailed logging.
 * 
 * @deprecated This function is deprecated and will be removed in future versions. Consider using the AdvanceParallel() function for parallel time step advancement.
 * 
 * @see AdvanceParallel(), IntegrateParallel()
 */
void LidDrivenCavity::Advance(bool verbose_advance)
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

    if(verbose_advance){
        std::cout << "Interior vorticity" << std::endl;
        Printmatrix(Nx,Ny,v);
        std::cout << "(Interior vorticity) Stream function" << std::endl;
        Printmatrix(Nx,Ny,s);
    }

    // Time advance vorticity
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            // v[IDX(i,j)] = v[IDX(i,j)] + dt*(
            v_new[IDX(i,j)] = v[IDX(i,j)] + dt*(
                ( (s[IDX(i+1,j)] - s[IDX(i-1,j)]) * 0.5 * dxi
                 *(v[IDX(i,j+1)] - v[IDX(i,j-1)]) * 0.5 * dyi)
              - ( (s[IDX(i,j+1)] - s[IDX(i,j-1)]) * 0.5 * dyi
                 *(v[IDX(i+1,j)] - v[IDX(i-1,j)]) * 0.5 * dxi)
              + nu * (v[IDX(i+1,j)] - 2.0 * v[IDX(i,j)] + v[IDX(i-1,j)])*dx2i
              + nu * (v[IDX(i,j+1)] - 2.0 * v[IDX(i,j)] + v[IDX(i,j-1)])*dy2i);
        }
    }

    if(verbose_advance){
        std::cout << "Next vorticity" << std::endl;
        Printmatrix(Nx,Ny,v);
        std::cout << "(Next voriticity) Stream function" << std::endl;
        Printmatrix(Nx,Ny,s);
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
    // cg->Solve(v, s, verbose);
    cg->Solve(v_new, s, verbose);


    if(verbose_advance){
        std::cout << "Laplace Operator (vorticity)" << std::endl;
        Printmatrix(Nx,Ny,v);
        std::cout << "Laplace Operator (stream function)" << std::endl;
        Printmatrix(Nx,Ny,s);
    }

}


/**
 * @brief Returns the grid spacing in the x-direction.
 * 
 * @return The grid spacing in the x-direction.
 */
double LidDrivenCavity::get_dx()
{
    return dx;
}


/**
 * @brief Returns the grid spacing in the y-direction.
 * 
 * @return The grid spacing in the y-direction.
 */
double LidDrivenCavity::get_dy()
{
    return dy;
}


/**
 * @brief Returns the total number of grid points in the domain.
 * 
 * @return The total number of grid points in the domain.
 */
int LidDrivenCavity::get_Npts()
{
    return Npts;
}


/**
 * @brief Returns the kinematic viscosity of the fluid.
 * 
 * @return The kinematic viscosity of the fluid.
 */
double LidDrivenCavity::get_nu()
{
    return nu;
}


/**
 * @brief Returns a pointer to the array storing the vorticity field.
 * 
 * @return A pointer to the array storing the vorticity field.
 */
double* LidDrivenCavity::get_v()
{
    return v;
}


/**
 * @brief Returns a pointer to the array storing the stream function field.
 * 
 * @return A pointer to the array storing the stream function field.
 */
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
        Advance(false);
    }
}