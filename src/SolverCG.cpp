#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <mpi.h>
#include <omp.h>
using namespace std;

#include <cblas.h>

#include "SolverCG.h"

/**
 * @brief Natural ordering based on global x and y coordinates
 */
#define IDX(I,J) ((J)*Nx + (I))
#define IDX_local(I,J) ((J)*Nx + (I)) // Nx_local same as Nx in SolverCG Class


/**
 * @brief Prints the content of a 2D matrix to the standard output.
 * @param nx The number of columns in the matrix.
 * @param ny The number of rows in the matrix.
 * @param A Pointer to the 1D array representing the matrix.
*/
void SolverCG::Printmatrix(int nx, int ny, double* A) {
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
 * @brief Initializes the SolverCG conjugate gradient solver solver for parallel execution.
 * 
 * Allocates memory for local residual vector (r), local conjugate direction vector (p), local preconditioned vector (z),
 * local temporary vector (t) as well as buffers required for inter-domain communication.
 * Initializes all arrays to prevent arbitrary and unpredictable values.
 * 
 * @param pNX The number of local grid points in x-direction
 * @param pNy The number of local grid points in y-direction
 * @param pdx Grid spacing in x-direction
 * @param pdy Grid spacing in y-direction
 * 
 * 
 */
SolverCG::SolverCG(int pNx, int pNy, double pdx, double pdy)
{
    dx = pdx;
    dy = pdy;
    Nx = pNx;
    Ny = pNy;
    int n = Nx*Ny;
    r = new double[n];
    p = new double[n];
    z = new double[n];
    t = new double[n]; //temp

    for(int i=0;i<n;i++){
        r[i] = 0;
        p[i] = 0;
        z[i] = 0;
        t[i] = 0;
    }


    // iter_max = 2;
    iter_max = 5000;

    // debug = true;
    debug = false;

    // Retrieve current rank and process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);

    buffer_up_send = new double[Nx];
    buffer_down_send = new double[Nx];
    buffer_up_recv = new double[Nx];
    buffer_down_recv = new double[Nx];
    for (int i=0;i<Nx;i++) {
        buffer_up_send[i] = 0;
        buffer_down_send[i] = 0;
        buffer_up_recv[i] = 0;
        buffer_down_recv[i] = 0;
    }    

    buffer_left_send = new double[Ny];
    buffer_right_send = new double[Ny];
    buffer_left_recv = new double[Ny];
    buffer_right_recv = new double[Ny];
    for (int j=0;j<Ny;j++) {
        buffer_left_send[j] = 0;
        buffer_right_send[j] = 0;
        buffer_left_recv[j] = 0;
        buffer_right_recv[j] = 0;
    }  

    dx2i = 1.0/dx/dx;
    dy2i = 1.0/dy/dy;
    factor2 = 1/(2.0*(dx2i + dy2i));

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
void SolverCG::SetNeighbour(int rank_up, int rank_down, int rank_left, int rank_right){
    this->rank_up = rank_up;
    this->rank_down = rank_down;
    this->rank_left = rank_left;
    this->rank_right = rank_right;

    i_start = 0;
    i_end = Nx;
    j_start = 0;
    j_end = Ny;
    if (this->rank_up != -2) {j_end--;}
    if (this->rank_down != -2) {j_start++;}
    if (this->rank_left != -2) {i_start++;}
    if (this->rank_right != -2) {i_end--;}

}

/**
 * @brief (Purely for debugging purposes) Configures the coordinates of the local domain relative to the global domain
 * and initialises memory at root process for assembling global domain
 * 
 * Global position can be computed with offset details, this faciliates gathering of all 
 * local domains for debugging purposes. Global matrices are initialised as zero
 * to prevent arbitrary and unpredictable values
 * 
 * @warning Will be removed future versions
 * 
 * @see GatherDomain()
*/
void SolverCG::SetOffset(int offset_x, int offset_y, int Nx_global, int Ny_global){
    this->offset_x = offset_x;
    this->offset_y = offset_y;
    this->Nx_global = Nx_global;
    this->Ny_global = Ny_global;   
    this->Npts = Nx_global*Ny_global;
}

/**
 * @brief Cleans up memory allocated for residue, conjugate direction, preconditioned, and temporary vector
*/
SolverCG::~SolverCG()
{
    delete[] r;
    delete[] p;
    delete[] z;
    delete[] t;

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
}

/**
 * @brief Solve a linear system using the Conjugate Gradient method in parallel.
 * 
 * 
 * @param b Pointer to the right-hand side vector (Vorticity in this code).
 * @param x Pointer to the solution vector (Stream function in this code).
 * @param verbose Whether to display convergence details during simulation
 * 
 * @details This method solves a linear system Ax = b using the Conjugate Gradient (CG) method
 * in parallel. It initializes necessary variables, computes error tolerance, and performs
 * iterative CG steps until convergence or reaching the maximum number of iterations.
 * 
 * @note Program terminates if convergence is not reached within the maximum number of iterations
 * 
 * @see ApplyOperatorParallel(), PreconditionParallel(), ImposeBCParallel()
*/
void SolverCG::SolveParallel(double* b, double* x, bool verbose) {
    unsigned int n = Nx*Ny;
    int k;
    double alpha_global, alpha1_global,alpha2_global;
    double beta_global, beta1_global, beta2_global;

    double eps;
    double tol = 0.001;

    eps = ComputeErrorGlobalParallel(n,b);
    if (eps < tol*tol) {
        std::fill(x, x+n, 0.0);
        cout << "Norm is " << eps << endl;
        return;
    }


    DomainInterComunnication(p); // Remove if unnecessary
    ApplyOperator(x, t);

    cblas_dcopy(n, b, 1, r, 1);        // r_0 = b (i.e. b)
    ImposeBCParallel(r);

    cblas_daxpy(n, -1.0, t, 1, r, 1);
    PreconditionParallel(r, z);
    DomainInterComunnication(z); // Remove if unnecessary


    cblas_dcopy(n, z, 1, p, 1);        // p_0 = r_0

    DomainInterComunnication(p); // Remove if unnecessary

    k = 0;
    do {
        k++;
        // Perform action of Nabla^2 * p
        ApplyOperator(p, t);
        DomainInterComunnication(t);

        alpha1_global = ComputeDotGlobalParallel(n,t,p);
        alpha2_global = ComputeDotGlobalParallel(n,r,z);
        alpha_global = alpha2_global/alpha1_global;

        // One single new search direction for all processes (beta_global)
        beta1_global  = alpha2_global;


        // cblas_daxpy(n,  alpha_global, p, 1, x, 1);  // x_{k+1} = x_k + alpha_k p_k
        // cblas_daxpy(n, -alpha_global, t, 1, r, 1); // r_{k+1} = r_k - alpha_k A p_k
        #pragma omp for
        for (unsigned int i = 0; i < n; ++i) {
            x[i] += alpha_global * p[i];
            r[i] -= alpha_global * t[i];
        }

        DomainInterComunnication(x); // Remove if unnecessary
        DomainInterComunnication(r); // Remove if unnecessary

        eps = ComputeErrorGlobalParallel(n,r);

        if (eps < tol*tol) {
            break;
        }
        PreconditionParallel(r, z);

        beta2_global = ComputeDotGlobalParallel(n,r,z);
        beta_global = beta2_global / beta1_global;

        // cblas_dcopy(n, z, 1, t, 1);
        // cblas_daxpy(n, beta_global, p, 1, t, 1);
        // cblas_dcopy(n, t, 1, p, 1);
        #pragma omp for 
        for (unsigned int i=0; i<n; i++) {
            t[i]=z[i];
            t[i]=z[i] + beta_global*p[i];
            p[i]=t[i];
        }

        DomainInterComunnication(p); // Remove if unnecessary
    } while (k < iter_max); // Set a maximum number of iterations


    if (rank==root && k == iter_max && debug==false) {
        cout << "FAILED TO CONVERGE" << endl;
        exit(-1);
    }

    if(rank==root && verbose){
        cout << "Converged in " << k << " iterations. eps = " << eps << endl;
    }
}


/**
 * @brief Apply the Laplace operator to the input vector.
 * 
 * @param in pointer to the input vector
 * @param out pointer to the operator vector
 * 
 * 
 * @details This method applies the Laplace operator to the input vector
 * representing the linear system, producing the output vector.
*/ 
void SolverCG::ApplyOperator(double* in, double* out) {
    // Assume ordered with y-direction fastest (column-by-column)
    // double dx2i = 1.0/dx/dx;
    // double dy2i = 1.0/dy/dy;
    // // int jm1 = 0, jp1 = 2;

    int index_local;
    #pragma omp for collapse(2) private(index_local)
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 1; i < Nx - 1; ++i) {
            index_local = IDX(i,j);

            // out[index_local] = ( -      in[index_local-1]
            //                     + 2.0*in[index_local]
            //                     -     in[index_local+1])*dx2i
            //                 + ( -     in[index_local-Nx]
            //                     + 2.0*in[index_local]
            //                     -     in[index_local+Nx])*dy2i;

            out[index_local] = (-in[index_local-1]-in[index_local+1])*dx2i + (2.0*dx2i+2.0*dy2i)*in[index_local] + (-in[index_local-Nx]-in[index_local+Nx])*dy2i;
            
        }
        // jm1++;
        // jp1++;
    }

}

/**
 * @brief Apply parallel preconditioning to the input vector.
 * 
 * @details This method applies parallel preconditioning to the input vector,
 * producing the output vector. Preconditioning is a technique used to improve
 * the convergence of iterative methods by transforming the linear system.
 * 
 * @param in Pointer to the input vector.
 * @param out Pointer to the output vector.
 */
void SolverCG::PreconditionParallel(double* in, double* out) {
    int i, j;
    // double dx2i = 1.0/dx/dx;
    // double dy2i = 1.0/dy/dy;
    // double factor = 2.0*(dx2i + dy2i);
    // double factor2 = 1/(2.0*(dx2i + dy2i));


    int i_start_temp = 1;
    int i_end_temp = Nx-1;
    int j_start_temp = 1;
    int j_end_temp = Ny-1;


    if (rank_up != -2) {j_end_temp++;}
    if (rank_down != -2) {;j_start_temp--;}
    if (rank_left != -2) {i_start_temp--;}
    if (rank_right != -2) {i_end_temp++;}



    for (j = j_start_temp; j < j_end_temp; ++j) {
        for (i = i_start_temp; i < i_end_temp; ++i) {
            // out[IDX(i,j)] = in[IDX(i,j)]/factor;
            out[IDX(i,j)] = in[IDX(i,j)]*factor2;
        }
    }
    // Boundaries
    for (i = 0; i < Nx; ++i) {
        if(rank_down==-2){
            out[IDX(i, 0)] = in[IDX(i,0)];
        }
        if(rank_up==-2){
            out[IDX(i, Ny-1)] = in[IDX(i, Ny-1)];
        }
    }

    for (j = 0; j < Ny; ++j) {
        if(rank_left==-2){
            out[IDX(0, j)] = in[IDX(0, j)];
        }
        if(rank_right==-2){
            out[IDX(Nx - 1, j)] = in[IDX(Nx - 1, j)];
        }
    }
}


/**
 * @brief Impose boundary conditions on the input vector in parallel.
 * 
 * @details This method imposes boundary conditions on the input vector in parallel,
 * modifying the values of boundary elements to satisfy specified boundary conditions.
 * 
 * @param inout Pointer to the modified input vector.
 */
void SolverCG::ImposeBCParallel(double* inout) {
    // Boundaries
    for (int i = 0; i < Nx; ++i) {
        if(rank_down==-2){
            inout[IDX(i, 0)] = 0.0;
        }
        if(rank_up==-2){
            inout[IDX(i, Ny-1)] = 0.0;
        }
    }

    for (int j = 0; j < Ny; ++j) {
        if(rank_left==-2){
            inout[IDX(0, j)] = 0.0;
        }
        if(rank_right==-2){
            inout[IDX(Nx - 1, j)] = 0.0;
        }
    }

}


/**
 * @brief Manages inter-domain communication for the conjugate gradient solver.
 * 
 * This function handles communication between neighboring processes in the parallel 
 * domain decomposition of the Lid-Driven Cavity simulation. It exchanges boundary 
 * information between adjacent domains to ensure consistency in the simulation results.
 * 
 * @param A_local Pointer to the local domain array containing the data to be communicated.
*/
void SolverCG::DomainInterComunnication(double* A_local){

    int Nx_local = Nx;
    int Ny_local = Ny;

    // Top and bottom communication
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
        MPI_Send(buffer_up_send, Nx_local, MPI_DOUBLE, rank_up, 0, MPI_COMM_WORLD);
    }
    if(rank_down!=-2){
        MPI_Send(buffer_down_send, Nx_local, MPI_DOUBLE, rank_down, 0, MPI_COMM_WORLD);
    }

    if(rank_down!=-2){
        MPI_Recv(buffer_down_recv, Nx_local, MPI_DOUBLE, rank_down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
    if(rank_up!=-2){
        MPI_Recv(buffer_up_recv, Nx_local, MPI_DOUBLE, rank_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    // Storing communicated knowledge
    for(int i=0;i<Nx_local;i++){
        if (rank_up != -2) {A_local[IDX_local(i,Ny_local-1)] = buffer_up_recv[i];}
        if (rank_down != -2) {A_local[IDX_local(i,0)] = buffer_down_recv[i];}
    }

    // Left and right communication
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
        MPI_Send(buffer_left_send, Ny_local, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD);
    }
    if(rank_right!=-2){
        MPI_Send(buffer_right_send, Ny_local, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD);
    }

    if(rank_right!=-2){
        MPI_Recv(buffer_right_recv, Ny_local, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(rank_left!=-2){
        MPI_Recv(buffer_left_recv, Ny_local, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Storing communicated knowledge
    for(int j=0;j<Ny_local;j++){
        if (rank_left != -2) {A_local[IDX_local(0,j)] = buffer_left_recv[j];}
        if (rank_right != -2) {A_local[IDX_local(Nx_local-1,j)] = buffer_right_recv[j];}
    }       
}


/**
 * @brief (Primarily for debugging purposes, not used in the code) Gathers the local domain data from each process and sends it to the root process to assemble the global domain.
 * The function ensures proper alignment and communication between processes in a parallel environment.
 * 
 * @param A_local Pointer to the local domain array.
 * @param A_global Pointer to the global domain array.
 * 
 * @details This function is responsible for gathering the local domain data from each process and sending it to the root process
 * to assemble the global domain. It takes into account the local offsets and grid sizes to ensure proper alignment and communication
 * between processes in a parallel environment
 * 
 * @deprecated This function is deprecated and will be removed in future versions.
 */
void SolverCG::GatherDomain(double* A_local, double* A_global){

    int Nx_local = Nx;
    int Ny_local = Ny;

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
                        index_global = (j + recv_buf[3])*Nx_global + (i + recv_buf[2]);
                        index_local = (j)*recv_buf[0] + (i);
                        A_global[index_global] = A_local_temp_recv[index_local];
                    }
                }
                delete[] A_local_temp_recv;
            }
            else{
                
                for(int i=i_start;i<i_end;i++){
                    for(int j=j_start;j<j_end;j++){
                        index_global = (j + offset_y)*Nx_global + (i + offset_x);
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

/**
 * @brief Compute the global error of a vector in parallel.
 * 
 * @details This method computes the global error of a vector in parallel by summing
 * the squared values of its elements across all MPI ranks and then taking the square root
 * of the sum.
 * 
 * @param n The total number of elements in the vector.
 * @param a Pointer to the local error vector.
 * 
 * @return The global reduced error.
 * 
 * @see ComputeDotGlobalParallel()
 */
double SolverCG::ComputeErrorGlobalParallel(unsigned int n, double *a){
    double error = 0; // squared without norm yet
    double error_global = 0;
    // error = cblas_ddot(n, a, 1, a, 1); // the square root process can only be done after global ddot error is gathered

    #pragma omp parallel for collapse(2) reduction(+:error)
    for(int j=j_start;j<j_end;j++){
        for(int i=i_start;i<i_end;i++){
            error += a[IDX_local(i,j)]*a[IDX_local(i,j)];
        }
    }

    MPI_Allreduce(&error, &error_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    error_global = sqrt(error_global);
    return error_global;
}

/**
 * @brief Compute the global dot prduct of two vectors in parallel.
 * 
 * @details This method computes the global dot product of two vector in parallel by summing
 * the values of its elements across all MPI ranks.
 * 
 * @param n The total number of elements in the vector.
 * @param a Pointer to the first vector.
 * @param b Pointer to the second vector.
 * 
 * 
 * @return The global dot product.
 * 
 * @see ComputeErrorGlobalParallel()
 */
double SolverCG::ComputeDotGlobalParallel(unsigned int n, double *a, double *b){
    double value = 0; // squared without norm yet
    double value_global = 0;

    int index_local;
    #pragma omp parallel for collapse(2) reduction(+:value) private(index_local)
    for(int j=j_start;j<j_end;j++){
        for(int i=i_start;i<i_end;i++){
            index_local = IDX_local(i,j);
            value += a[index_local]*b[index_local];
        }
    }

    MPI_Allreduce(&value, &value_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return value_global;
}



























//////////////////////////////////////////    Legacy baseline code    //////////////////////////////////////////


/**
 * @brief Solve a linear system using the Conjugate Gradient method in serial.
 * 
 * 
 * @param b Pointer to the right-hand side vector (Vorticity in this code).
 * @param x Pointer to the solution vector (Stream function in this code).
 * @param verbose Whether to display convergence details during simulation
 * 
 * @details This method solves a linear system Ax = b using the Conjugate Gradient (CG) method
 * in serial. It initializes necessary variables, computes error tolerance, and performs
 * iterative CG steps until convergence or reaching the maximum number of iterations.
 * 
 * @deprecated This function is deprecated and will be removed in future versions. Consider using the SolveParallel() function for parallel cg solve.
 * 
 * @note Program terminates if convergence is not reached within the maximum number of iterations
 * 
 * @see ApplyOperatorParallel(), PreconditionParallel(), ImposeBCParallel()
*/
void SolverCG::Solve(double* b, double* x, bool verbose) {
    unsigned int n = Nx*Ny;
    int k;
    double alpha;
    double beta;
    double eps;
    double tol = 0.001;

    eps = cblas_dnrm2(n, b, 1);
    if (eps < tol*tol) {
        std::fill(x, x+n, 0.0);
        cout << "Norm is " << eps << endl;
        return;
    }
    ApplyOperator(x, t);
    cblas_dcopy(n, b, 1, r, 1);        // r_0 = b (i.e. b)
    ImposeBC(r);


    cblas_daxpy(n, -1.0, t, 1, r, 1);
    Precondition(r, z);

    cblas_dcopy(n, z, 1, p, 1);        // p_0 = r_0

    k = 0;
    do {
        k++;
        // Perform action of Nabla^2 * p
        ApplyOperator(p, t);

        alpha = cblas_ddot(n, t, 1, p, 1);  // alpha = p_k^T A p_k

        alpha = cblas_ddot(n, r, 1, z, 1) / alpha; // compute alpha_k
        beta  = cblas_ddot(n, r, 1, z, 1);  // z_k^T r_k

        cblas_daxpy(n,  alpha, p, 1, x, 1);  // x_{k+1} = x_k + alpha_k p_k
        cblas_daxpy(n, -alpha, t, 1, r, 1); // r_{k+1} = r_k - alpha_k A p_k

        eps = cblas_dnrm2(n, r, 1);

        if (eps < tol*tol) {
            break;
        }
        Precondition(r, z);

        beta = cblas_ddot(n, r, 1, z, 1) / beta;

        cblas_dcopy(n, z, 1, t, 1);
        cblas_daxpy(n, beta, p, 1, t, 1);
        cblas_dcopy(n, t, 1, p, 1);

    } while (k < iter_max); // Set a maximum number of iterations

    if (k == iter_max && debug==false) {
        cout << "FAILED TO CONVERGE" << endl;
        exit(-1);
    }

    if(verbose){
        cout << "Converged in " << k << " iterations. eps = " << eps << endl;
    }
}


/**
 * @brief Apply serial preconditioning to the input vector.
 * 
 * @details This method applies preconditioning to the input vector,
 * producing the output vector. Preconditioning is a technique used to improve
 * the convergence of iterative methods by transforming the linear system.
 * 
 * @param in Pointer to the input vector.
 * @param out Pointer to the output vector.
 * 
 * @deprecated This function is deprecated and will be removed in future versions. Consider using the PreconditionParallel() function for parallel preconditioning.
 * 
 */
void SolverCG::Precondition(double* in, double* out) {
    int i, j;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    double factor = 2.0*(dx2i + dy2i);
    for (i = 1; i < Nx - 1; ++i) {
        for (j = 1; j < Ny - 1; ++j) {
            out[IDX(i,j)] = in[IDX(i,j)]/factor;
        }
    }
    // Boundaries
    for (i = 0; i < Nx; ++i) {
        out[IDX(i, 0)] = in[IDX(i,0)];
        out[IDX(i, Ny-1)] = in[IDX(i, Ny-1)];
    }

    for (j = 0; j < Ny; ++j) {
        out[IDX(0, j)] = in[IDX(0, j)];
        out[IDX(Nx - 1, j)] = in[IDX(Nx - 1, j)];
    }
}


/**
 * @brief Impose boundary conditions on the input vector in serial.
 * 
 * @details This method imposes boundary conditions on the input vector in serial.
 * modifying the values of boundary elements to satisfy specified boundary conditions.
 * 
 * @param inout Pointer to the modified input vector.
 * 
 * @deprecated This function is deprecated and will be removed in future versions. Consider using the ImposeBCParallel() function for parallel boundary imposing
 */
void SolverCG::ImposeBC(double* inout) {
        // Boundaries
    for (int i = 0; i < Nx; ++i) {
        inout[IDX(i, 0)] = 0.0;
        inout[IDX(i, Ny-1)] = 0.0;
    }

    for (int j = 0; j < Ny; ++j) {
        inout[IDX(0, j)] = 0.0;
        inout[IDX(Nx - 1, j)] = 0.0;
    }

}
