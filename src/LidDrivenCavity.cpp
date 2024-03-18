#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <mpi.h>
#include <omp.h>

using namespace std;

#include <cblas.h>

#define IDX(I,J) ((J)*Nx + (I))
#define IDX_local(I_local,J_local) ((J_local)*Nx_local + (I_local))
#define IDX_local_buffer(I_local,J_local) ((J_local+1)*(Nx_local+2) + (I_local+1))



#include "LidDrivenCavity.h"
#include "SolverCG.h"

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

void LidDrivenCavity::PrintMatrix(int nsv, double* A) {
    cout.precision(4);
    for (int i = 0; i < nsv; ++i) {
        for (int j = 0; j < nsv; ++j) {
            cout << setw(13) << A[j*nsv+i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int LidDrivenCavity::Local2Global(int i_local, int j_local){
    int index_global = (j_local + offset_y + 1)*Nx + (i_local + offset_x + 1);
    return index_global;
}


// Function to check if a number is at the boundary of a grid
bool LidDrivenCavity::CheckBoundary(int i_local, int j_local)
{
    int index = Local2Global(i_local,j_local);
    int rows = Ny;
    int cols = Nx;
    return index < cols ||                // Top row
           index >= (rows - 1) * cols ||  // Bottom row
           index % cols == 0 ||           // Leftmost column
           index % cols == cols - 1;      // Rightmost column
}

void LidDrivenCavity::ScatterDomain(double* A_local, double* A_global)
{



    if(rank==root){
        for (int src = 1; src < Nprocs; ++src) {
            MPI_Send(A_global, Npts, MPI_DOUBLE, src, 0, MPI_COMM_WORLD);
        }
        int index_global, index_local;
        for (int i=0;i<Nx_local;i++){
            for (int j=0;j<Ny_local;j++){
                index_global = (j + 0)*Nx + (i + 0);
                index_local = IDX_local(i,j);                
                A_local[index_local] = A_global[index_global];
            }
        }
    } 
    else{
        double* A_global_temp = new double[Npts];
        int index_global, index_local;
        MPI_Recv(A_global_temp, Npts, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int offset_y_temp = offset_y;
        int offset_x_temp = offset_x;
        if(rank_down!=-2){offset_y_temp--;}
        if(rank_left!=-2){offset_x_temp--;}


        for (int i=0;i<Nx_local;i++){
            for (int j=0;j<Ny_local;j++){
                index_global = (j + offset_y_temp)*Nx + (i + offset_x_temp);
                index_local = IDX_local(i,j);                
                A_local[index_local] = A_global_temp[index_global];
            }
        }
        delete[] A_global_temp;
    }

}


void LidDrivenCavity::GatherDomain(double* A_local, double* A_global){

    int Nx_local_temp = Nx_local;
    int Ny_local_temp = Ny_local;
    int i_start = 0;
    int i_end = Nx_local;
    int j_start = 0;
    int j_end = Ny_local;
    if (rank_up != -2) {Ny_local_temp--;j_end--;}
    if (rank_down != -2) {Ny_local_temp--;j_start++;}
    if (rank_left != -2) {Nx_local_temp--;i_start++;}
    if (rank_right != -2) {Nx_local_temp--;i_end--;}

    int Npts_local_temp = Nx_local_temp*Ny_local_temp;
    double* A_local_temp = new double[Npts_local_temp];
    for(int i=0;i<Npts_local_temp;i++){A_local_temp[i]=0;}

    int k = 0;
    for(int j=j_start;j<j_end;j++){
        for(int i=i_start;i<i_end;i++){
            A_local_temp[k] = A_local[IDX_local(i,j)];
            k++;
        }
    }

    if(rank!=root){

        int send_buf[4] = {Nx_local_temp, Ny_local_temp, offset_x, offset_y};
        MPI_Send(send_buf,4,MPI_INT,root,0,MPI_COMM_WORLD);
        // Send A_local_temp to the root
        MPI_Send(A_local_temp, Npts_local_temp, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
    }
    else{
        int index_global, index_local;
        for (int src = 0; src < Nprocs; ++src) {

            if (src != root) {
                int recv_buf[4];
                MPI_Recv(recv_buf, 4, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                double* A_local_temp_recv = new double[recv_buf[0]*recv_buf[1]];

                MPI_Recv(A_local_temp_recv, recv_buf[0]*recv_buf[1], MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int i=0;i<recv_buf[0];i++){
                    for(int j=0;j<recv_buf[1];j++){
                        index_global = (j + recv_buf[3])*Nx + (i + recv_buf[2]);
                        index_local = (j)*recv_buf[0] + (i);
                        A_global[index_global] = A_local_temp_recv[index_local];

                    }
                }
                delete[] A_local_temp_recv;
            }
            else{
                
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

void LidDrivenCavity::SetThreads(int Nt)
{
    this->Nt = Nt;
}

void LidDrivenCavity::SetNeighbour(int rank_up, int rank_down, int rank_left, int rank_right)
{
    this->rank_up = rank_up;
    this->rank_down = rank_down;
    this->rank_left = rank_left;
    this->rank_right = rank_right;

}

void LidDrivenCavity::DomainDecomposition()
{

    // Retrieve current rank and process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);


    parallel = Nprocs > 1;
    // if(parallel){
        Nprocs_sqrt = sqrt(Nprocs);
        Nx_remainder = Nx % Nprocs_sqrt;
        Ny_remainder = Ny % Nprocs_sqrt;

        if(rank % Nprocs_sqrt < Nx_remainder){
            Nx_local = Nx / Nprocs_sqrt + 1;
            offset_x = Nx_local * (rank % Nprocs_sqrt);
        }
        else{
            Nx_local = Nx / Nprocs_sqrt;
            offset_x = Nx_local * (rank % Nprocs_sqrt) + Nx_remainder;
        }

        if(rank / Nprocs_sqrt < Ny_remainder){
            Ny_local = Ny / Nprocs_sqrt + 1;
            offset_y = Ny_local * (rank / Nprocs_sqrt);
        }
        else{
            Ny_local = Ny / Nprocs_sqrt;
            offset_y = Ny_local * (rank / Nprocs_sqrt) + Ny_remainder ;
        }

        // Account for shared rows/columns of neighbouring ranks
        if (rank_up != -2) {Ny_local++;}
        if (rank_down != -2) {Ny_local++;}
        if (rank_left != -2) {Nx_local++;}
        if (rank_right != -2) {Nx_local++;}

        Npts_local = Nx_local*Ny_local;

        // std::cout << "Rank: " << rank  << " | Nx_local: " << Nx_local << " | offset_x: " << offset_x << " | Ny_local: " << Ny_local << " | offset_y: " << offset_y << std::endl;

}

void LidDrivenCavity::InitialiseParallel()
{
    v_local = new double[Npts_local];
    s_local = new double[Npts_local];
    v_next_local = new double[Npts_local];
    cg = new SolverCG(Nx_local,Ny_local,dx,dy);
    cg->SetNeighbour(rank_up,rank_down,rank_left,rank_right);
    cg->SetOffset(offset_x,offset_y,Nx,Ny);

    // Initialize elements for safety
    for (int i=0;i<Npts_local;i++) {
        v_local[i] = 0;
        s_local[i] = 0;
        v_next_local[i] = 0;
    }    

    if(rank==root){
        u0 = new double[Npts];
        u1 = new double[Npts];
        for (int i=0;i<Npts;i++) {u0[i] = 0;u1[i] = 0;}
    }
    if(rank==root){
        v = new double[Npts];
        s = new double[Npts];
        for (int i=0;i<Npts;i++) {v[i] = 0;s[i]=0;}
    }

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

void LidDrivenCavity::CreateU(){
    u0 = new double[Npts];
    u1 = new double[Npts];
    for (int i=0;i<Npts;i++) {u0[i] = 0;u1[i] = 0;}

}

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
        if(t==-1){verbose_advance=true;}
        else{verbose_advance=false;}
        AdvanceParallel(verbose_advance);
    }
    GatherDomain(v_local,v);
    GatherDomain(s_local,s);
}

void LidDrivenCavity::AdvanceParallel(bool verbose_advance){

    ComputeBoundaryVorticityParallel();
    ComputeInteriorVorticityParallel();

    // GatherDomain(v_local,v);
    if(rank==root&&verbose_advance){
        std::cout << "(Interior vorticity) vorticity matrix" << std::endl;
        Printmatrix(Nx, Ny, v);
    }

    DomainInterComunnication(v_local);
    ComputeNextVorticityParallel();
    DomainInterComunnication(v_local);  

    // GatherDomain(v_local,v);
    // GatherDomain(s_local,s);
    if(rank==root&&verbose_advance){
        std::cout << "(Next vorticity) vorticity matrix" << std::endl;
        Printmatrix(Nx, Ny, v);
        std::cout << "(Next vorticity) stream function matrix" << std::endl;
        Printmatrix(Nx, Ny, s);
    } 

    ComputeLaplaceOperatorParallel();

    // GatherDomain(v_local,v);
    // GatherDomain(s_local,s);
    if(rank==root&&verbose_advance){
        std::cout << "(Laplace operator) vorticity matrix" << std::endl;
        Printmatrix(Nx, Ny, v);
        std::cout << "(Laplace operator) stream function matrix" << std::endl;
        Printmatrix(Nx, Ny, s);
    }   

}


void LidDrivenCavity::ComputeBoundaryVorticityParallel()
{

    double dyi  = 1.0/dy;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;

    // Top
    if (rank_up==-2){
        for (int i=1;i<Nx_local-1;++i){
            v_local[IDX_local(i,Ny_local-1)] = 2.0 * dy2i * (s_local[IDX_local(i,Ny_local-1)] - s_local[IDX_local(i,Ny_local-2)])- 2.0 * dyi*U;
        }
    }

    // Bottom
    if (rank_down==-2){
        for (int i=1;i<Nx_local-1;++i){
            v_local[IDX_local(i,0)] = 2.0 * dy2i * (s_local[IDX_local(i,0)] - s_local[IDX_local(i,1)]);
        }
    }

    // Left
    if (rank_left==-2){
        for (int j=1;j<Ny_local-1;++j){
            v_local[IDX_local(0,j)]    = 2.0 * dx2i * (s_local[IDX_local(0,j)] - s_local[IDX_local(1,j)]);
        }
    }

    // Right
    if (rank_right==-2){
        for (int j=1;j<Ny_local-1;++j){
            v_local[IDX_local(Nx_local-1,j)] = 2.0 * dx2i * (s_local[IDX_local(Nx_local-1,j)]    - s_local[IDX_local(Nx_local-2,j)]);
        }
    }


    // std::cout << "Compute Boundary" << std::endl;
}

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

void LidDrivenCavity::ComputeNextVorticityParallel(){
    double dxi  = 1.0/dx;
    double dyi  = 1.0/dy;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    // Time advance vorticity
    for (int i=1;i<Nx_local-1;i++){
        for (int j=1;j<Ny_local-1;j++){ // Check whether to use another variable or not
            v_local[IDX_local(i,j)] = v_local[IDX_local(i,j)] + dt*(
                ( (s_local[IDX_local(i+1,j)] - s_local[IDX_local(i-1,j)]) * 0.5 * dxi
                 *(v_local[IDX_local(i,j+1)] - v_local[IDX_local(i,j-1)]) * 0.5 * dyi)
              - ( (s_local[IDX_local(i,j+1)] - s_local[IDX_local(i,j-1)]) * 0.5 * dyi
                 *(v_local[IDX_local(i+1,j)] - v_local[IDX_local(i-1,j)]) * 0.5 * dxi)
              + nu * (v_local[IDX_local(i+1,j)] - 2.0 * v_local[IDX_local(i,j)] + v_local[IDX_local(i-1,j)])*dx2i
              + nu * (v_local[IDX_local(i,j+1)] - 2.0 * v_local[IDX_local(i,j)] + v_local[IDX_local(i,j-1)])*dy2i);
        }
    }
}


void LidDrivenCavity::ComputeLaplaceOperatorParallel()
{
    // Solve Poisson problem
    cg->SolveParallel(v_local, s_local, verbose);
 
}

void LidDrivenCavity::DomainInterComunnication(double* A_local){

    // Top and bottom communication
    // double buffer_up_send[Nx_local], buffer_down_send[Nx_local];
    // double buffer_up_recv[Nx_local], buffer_down_recv[Nx_local];

    if(rank_up!=-2){
        int j_second_upper = (Ny_local-1)-1;
        for(int i=0;i<Nx_local;i++){
            buffer_up_send[i] = A_local[IDX_local(i,j_second_upper)];
        }
    }
    if(rank_down!=-2){
        int j_second_lower = (0)+1;
        for(int i=0;i<Nx_local;i++){
            buffer_down_send[i] = A_local[IDX_local(i,j_second_lower)];
        }
    }

    // Send and receive
    if(rank_up!=-2){
        // MPI_Send(&buffer_up_send, Nx_local, MPI_DOUBLE, rank_up, 0, MPI_COMM_WORLD);
        MPI_Send(buffer_up_send, Nx_local, MPI_DOUBLE, rank_up, 0, MPI_COMM_WORLD);
    }
    if(rank_down!=-2){
        // MPI_Send(&buffer_down_send, Nx_local, MPI_DOUBLE, rank_down, 0, MPI_COMM_WORLD);
        MPI_Send(buffer_down_send, Nx_local, MPI_DOUBLE, rank_down, 0, MPI_COMM_WORLD);
    }

    if(rank_down!=-2){
        // MPI_Recv(&buffer_down_recv, Nx_local, MPI_DOUBLE, rank_down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(buffer_down_recv, Nx_local, MPI_DOUBLE, rank_down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(rank_up!=-2){
        // MPI_Recv(&buffer_up_recv, Nx_local, MPI_DOUBLE, rank_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(buffer_up_recv, Nx_local, MPI_DOUBLE, rank_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    // Storing communicated knowledge
    for(int i=0;i<Nx_local;i++){
        if (rank_up != -2) {A_local[IDX_local(i,Ny_local-1)] = buffer_up_recv[i];}
        if (rank_down != -2) {A_local[IDX_local(i,0)] = buffer_down_recv[i];}
    }

    // Left and right communication
    // Allocate memory for communication buffers
    // double buffer_left_send[Ny_local], buffer_right_send[Ny_local];
    // double buffer_left_recv[Ny_local], buffer_right_recv[Ny_local];

    if(rank_left!=-2){
        int i_second_left = (0) + 1;
        for(int j=0;j<Ny_local;j++){
            buffer_left_send[j] = A_local[IDX_local(i_second_left,j)];
        }
    }
    if(rank_right!=-2){
        int i_second_right = (Nx_local-1) - 1;
        for(int j=0;j<Ny_local;j++){
            buffer_right_send[j] = A_local[IDX_local(i_second_right,j)];
        }
    }

    // Send and receive
    if(rank_left!=-2){
        // MPI_Send(&buffer_left_send, Ny_local, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD);
        MPI_Send(buffer_left_send, Ny_local, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD);
    }
    if(rank_right!=-2){
        // MPI_Send(&buffer_right_send, Ny_local, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD);
        MPI_Send(buffer_right_send, Ny_local, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD);
    }

    if(rank_right!=-2){
        // MPI_Recv(&buffer_right_recv, Ny_local, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(buffer_right_recv, Ny_local, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(rank_left!=-2){
        // MPI_Recv(&buffer_left_recv, Ny_local, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(buffer_left_recv, Ny_local, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Storing communicated knowledge
    for(int j=0;j<Ny_local;j++){
        if (rank_left != -2) {A_local[IDX_local(0,j)] = buffer_left_recv[j];}
        if (rank_right != -2) {A_local[IDX_local(Nx_local-1,j)] = buffer_right_recv[j];}
    }       
}








































void LidDrivenCavity::Initialise()
{
    CleanUp();

    v   = new double[Npts]();
    s   = new double[Npts]();
    // tmp = new double[Npts]();
    cg  = new SolverCG(Nx, Ny, dx, dy);
    CreateU();
}

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

void LidDrivenCavity::WriteSolution(std::string file)
{

    for(int i=0;i<Npts;i++){u0[i] = 0; u1[i] = 0;}
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



    // delete[] u0;
    // delete[] u1;
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
    cout << "Number of threads: " << Nt << endl;
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
        // delete[] tmp;
        delete cg;
    }
}


void LidDrivenCavity::UpdateDxDy()
{
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    Npts = Nx * Ny;
}


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
        PrintMatrix(Nx,v);
        std::cout << "(Interior vorticity) Stream function" << std::endl;
        PrintMatrix(Nx,s);
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
    cg->Solve(v, s, verbose);

    if(verbose_advance){
        std::cout << "Laplace Operator (vorticity)" << std::endl;
        PrintMatrix(Nx,v);
        std::cout << "Laplace Operator (stream function)" << std::endl;
        PrintMatrix(Nx,s);
    }

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
        Advance(false);
    }
}