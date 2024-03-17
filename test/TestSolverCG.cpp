/*
Test cases for LidDrivenCavity class

Solver (Serial) Check
1. Arbitrary matrix check: Vorticity and stream function as arbitrary matrix
2. Real case check: real and realistic v and s taken from LidDrivenCavity class

Solver (Parallel) Check

*/

#include "../src/SolverCG.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/test/unit_test.hpp>
#include <mpi.h>
#include <fstream>
BOOST_AUTO_TEST_SUITE(solver_cg)

void PrintMatrix(double* A, int nx, int ny) {
    std::cout.precision(4);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            std::cout << std::setw(10) << A[j*nx+i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Function to load data from a file into a vector
void LoadDataFromFile(const std::string& file_path, double* A, int Npts) {
    std::ifstream infile(file_path);
    
    // Check if the file is opened successfully
    if (!infile) {
        std::cerr << "Error: Unable to open the file '" << file_path << "' for reading!" << std::endl;
        return;
    }

    // Load data from the file
    for (int i = 0; i < Npts; i++) {
        infile >> A[i];
    }

    // Close the file
    infile.close();
}



// Test case for arbitrary vorticity and stream function
BOOST_AUTO_TEST_CASE(SolveArbitraryCheck)
{
    // Test cases parameters
    int Nx_test = 29;
    int Ny_test = 25;

    int Npts = Nx_test*Ny_test;
    double dx_test = 0.008;
    double dy_test = 0.01;

    double* v = new double[Npts];
    double* s = new double[Npts];
    double* s_true = new double[Npts];
    bool verbose = false;

    // Create SolverCg environment
    SolverCG* solver_cg  = new SolverCG(Nx_test,Ny_test,dx_test,dy_test);

    for(int i=0;i<Npts;i++){
        s[i] = i*2e-5-5e-5;
        v[i] = (i % Nx_test);
    }

    // PrintMatrix(s,Nx_test,Ny_test);

    solver_cg->Solve(v,s,verbose);
    // PrintMatrix(s,Nx_test,Ny_test);

    // Open the file for reading
    std::string file_path = "test/data/arbitrary1.txt";
    LoadDataFromFile(file_path,s_true,Npts);

    // Check if solves correctly
    double tolerance = 1e-6;
    bool all_within = true;
    for (int i=0;i<Npts;i++){
        if(std::abs(s[i]-s_true[i])>tolerance){
            all_within = false;
        }
    }
    BOOST_CHECK_MESSAGE(all_within, "Solved stream function not within tolerance");
    BOOST_TEST_MESSAGE("\n  Test results: Conjugate gradient solver solves correctly\n");  // Print a message indicating the test case is successful

    delete solver_cg;
    delete[] v;
    delete[] s;
    delete[] s_true;
}


// Test case for real vorticity and stream function
BOOST_AUTO_TEST_CASE(SolveRealCheck)
{

    double Lx_test = 2.5;
    double Ly_test = 1.3;
    // Grid points for test case
    int Nx_test = 39;
    int Ny_test = 47;
    int Npts = Nx_test*Ny_test;

    double dx_test = Lx_test/(Nx_test+1);
    double dy_test = Ly_test/(Ny_test+1);
    double* v = new double[Npts];
    double* s = new double[Npts];
    double* v_true = new double[Npts];
    double* s_true = new double[Npts];
    bool verbose = false;

    // Open the file for reading
    std::string file_path = "test/data/realv_input1.txt";
    LoadDataFromFile(file_path,v,Npts);
    file_path = "test/data/reals_input1.txt";
    LoadDataFromFile(file_path,s,Npts);
    file_path = "test/data/realv_output1.txt";
    LoadDataFromFile(file_path,v_true,Npts);
    file_path = "test/data/reals_output1.txt";
    LoadDataFromFile(file_path,s_true,Npts);


    // Create SolverCg environment
    SolverCG* solver_cg  = new SolverCG(Nx_test,Ny_test,dx_test,dy_test);
    solver_cg->Solve(v,s,verbose);

    // Check if solves correctly
    double tolerance = 1e-6;
    bool all_within = true;
    for (int i=0;i<Npts;i++){
        if(std::abs(s[i]-s_true[i])>tolerance){
            all_within = false;
        }
    }
    BOOST_CHECK_MESSAGE(all_within, "Solved stream function not within tolerance");

}

/*
void ScatterDomain(double* A_local, double* A_global, int Nx_local, int Ny_local,
                   int offset_x, int offset_y)
{

    int root = 0;

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

void GatherDomain(double* A_local, double* A_global, int Nx_local, int Ny_local, int Nx, int Ny,
                  int Nprocs, int rank, int rank_up, int rank_down, int rank_left, int rank_right,
                  int offset_x, int offset_y){

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

    int IDX_local;

    int k = 0;
    for(int j=j_start;j<j_end;j++){
        for(int i=i_start;i<i_end;i++){
            IDX_local = (j)*Nx_local + (i);
            A_local_temp[k] = A_local[IDX_local];
            k++;
        }
    }

    int root = 0;
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
                        IDX_local = (j)*Nx_local + (i);
                        A_global[index_global] = A_local[IDX_local];
                    
                    }
                }
            }
                
        }
    }
    // Synchronize all processes
    MPI_Barrier(MPI_COMM_WORLD);

    delete[] A_local_temp;

}

// Test case for arbitrary vorticity and stream function
BOOST_AUTO_TEST_CASE(SolveArbitraryParallelCheck)
{
    // Test cases parameters
    int Nx = 29;
    int Ny = 25;

    int Npts = Nx*Ny;
    double dx = 0.008;
    double dy = 0.01;

    double* v = new double[Npts];
    double* s = new double[Npts];
    double* s_true = new double[Npts];
    bool verbose = false;

    for(int i=0;i<Npts;i++){
        s[i] = i*2e-5-5e-5;
        v[i] = (i % Nx);
    }

    int rank, Nprocs;
    int Nprocs_sqrt, Nx_remainder, Ny_remainder;
    int Nx_local, Ny_local, offset_x, offset_y;
    int Npts_local;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);

    // Create a cartesian grid of n x n to compute
    MPI_Comm domain_local;
    int dims[2] = {0, 0};
    int periods[2] = {0, 0};
    int coords[2];
    int reorder = 0; // Do not allow reordering for now

    MPI_Dims_create(Nprocs, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &domain_local); // Create the Cartesian communicator
    MPI_Cart_coords(domain_local, rank, 2, coords); // Get the coordinates of the current process in the Cartesian grid

    int rank_left, rank_right, rank_up, rank_down; // Ranks of neighboring processes
    // Get the ranks of the neighboring processes
    MPI_Cart_shift(domain_local, 0, 1, &rank_down, &rank_up);
    MPI_Cart_shift(domain_local, 1, 1, &rank_left, &rank_right);

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





    // Create SolverCg environment
    SolverCG* solver_cg  = new SolverCG(Nx_local,Ny_local,dx,dy);

    

    // Open the file for reading
    std::string file_path = "test/data/arbitrary1.txt";
    LoadDataFromFile(file_path,s_true,Npts);

    // Check if solves correctly
    double tolerance = 1e-6;
    bool all_within = true;
    for (int i=0;i<Npts;i++){
        if(std::abs(s[i]-s_true[i])>tolerance){
            all_within = false;
        }
    }
    BOOST_CHECK_MESSAGE(all_within, "Solved stream function not within tolerance");
    BOOST_TEST_MESSAGE("\n  Test results: Conjugate gradient solver solves correctly\n");  // Print a message indicating the test case is successful

    delete solver_cg;
    delete[] v;
    delete[] s;
    delete[] s_true;
}
*/


BOOST_AUTO_TEST_SUITE_END()

