/*
Test cases for LidDrivenCavity class

Program arguments
1. DomainCheck: Check correct mesh size
2. ReynoldsCheck: Check correct Reynolds number and Nusselt number specified

Numerical simulation
3. InitialConditionCheck: Check if initial condition (v,s,u0,u1) is applied correctly
4. BoundaryConditionCheck: Check if boundary condition is been applied correctly at the start, middle, and end of simulation
*/


#include "../src/LidDrivenCavity.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>
#include <memory>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(MyTests)

// Test case for correct mesh size specified
BOOST_AUTO_TEST_CASE(DomainCheck)
{
    // Domain length for test case
    double Lx_test = 2.0;
    double Ly_test = 1.5;
    // Grid points for test case
    int Nx_test = 39;
    int Ny_test = 121;

    // Create LidDrivenCavity environment
    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetDomainSize(Lx_test,Ly_test);
    solver->SetGridSize(Nx_test,Ny_test);

    // Check if mesh size is equivalent
    BOOST_CHECK_MESSAGE(solver->get_dx() == Lx_test / (Nx_test - 1), "Mesh size in x-direction does not match");
    BOOST_CHECK_MESSAGE(solver->get_dy() == Ly_test / (Ny_test - 1), "Mesh size in y-direction does not match");
    BOOST_TEST_MESSAGE("\n  Test results: Mesh size matches correctly\n");  // Print a message indicating the test case is successful

    delete solver;
}

// Test case for proper Reynolds and Nusselt number specified
BOOST_AUTO_TEST_CASE(ReynoldsCheck)
{
    // Reynolds and Nusselt number for test case
    double Re_test = 79.0;
    double nu_test = 1.0/Re_test;

    // Create LidDrivenCavity environment
    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetReynoldsNumber(Re_test);

    // Check if Nu from solver and Nu_test is equivalent
    BOOST_CHECK_MESSAGE(solver->get_nu() == nu_test, "Nusselt number does not match");
    BOOST_TEST_MESSAGE("\n  Test results: Reynolds and Nusselt number specified correctly\n");  // Print a message indicating the test case is successful

    delete solver;
}

// Check if all elements are zero
bool CheckZeroMatrix(const double* A, int size) {
    bool all_zero;
    all_zero = boost::algorithm::all_of(A, A + size, [](double val) { return val == 0.0; });
    return all_zero;
}

// Extract the inner rectangle (Nx-1)x(Ny-1) sub matrix B from matrix A
void ExtractInnerRectangle(double* A, double* B, int Nx, int Ny) {
    for(int i=1;i<Nx-1;i++){
        for(int j=1;j<Ny-1;j++){
            B[(j-1)*(Nx-2)+(i-1)] = A[j*Nx+i];
        }
    }
}

// Test case for initial condition (fluid at t=0 is at rest)
BOOST_AUTO_TEST_CASE(InitialConditionCheck)
{
    // Domain length for test case
    double Lx_test = 2.0;
    double Ly_test = 1.7;
    // Grid points for test case
    int Nx_test = 57;
    int Ny_test = 151;

    // Create LidDrivenCavity environment
    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetDomainSize(Lx_test,Ly_test);
    solver->SetGridSize(Nx_test,Ny_test);
    solver->Initialise();

    // Vorticity and stream function
    int Npts = solver->get_Npts();
    double* v = solver->get_v();
    double* s = solver->get_s();

    // Check if vorticity and stream function are zero
    BOOST_CHECK_MESSAGE(CheckZeroMatrix(v,Npts), "Vorticity matrix is not zero, fluid is not at rest");
    BOOST_CHECK_MESSAGE(CheckZeroMatrix(s,Npts), "Stream function matrix is not zero, fluid is not at rest"); 

    // velocity in x (u0) and y (u1) direction
    int Npts_inner = (Nx_test-2)*(Ny_test-2);
    double* u0 = solver->get_u0();
    double* u1 = solver->get_u1();
    double* u0_inner = new double[Npts_inner];
    double* u1_inner = new double[Npts_inner];
    // Exclude boundary condition elements from u0 and u1
    ExtractInnerRectangle(u0,u0_inner,Nx_test,Ny_test);
    ExtractInnerRectangle(u1,u1_inner,Nx_test,Ny_test);

    // // Check if inner velocity (excluding boundary condition)are zero
    BOOST_CHECK_MESSAGE(CheckZeroMatrix(u0_inner,Npts_inner), "U velocity matrix is not zero, fluid is not at rest");
    BOOST_CHECK_MESSAGE(CheckZeroMatrix(u1_inner,Npts_inner), "V velocity matrix is not zero, fluid is not at rest"); 

    BOOST_TEST_MESSAGE("\n  Test results: Initial condition applied correctly\n");  // Print a message indicating the test case is successful

    delete[] u0;
    delete[] u0_inner;
    delete[] u1;
    delete[] u1_inner;
    delete solver;
}

void ExtractRow(double* A, double* B, int Nx, int row){
    for(int i=0;i<Nx;i++){
        B[i] = A[row*Nx+i];
    }
}

void ExtractColumn(double* A, double* B, int Nx, int Ny, int col){
    for(int j=0;j<Ny;j++){
        B[j] = A[j*Nx+col];
    }
}

// Check if all elements are zero
bool CheckOneMatrix(const double* A, int size) {
    bool all_one;
    all_one = boost::algorithm::all_of(A, A + size, [](double val) { return val == 1.0; });
    return all_one;
}

// Check if all elements are zero
bool CheckZeroOneMatrix(const double* A, int size) {
    bool all_zero, last_one;
    all_zero = boost::algorithm::all_of(A, A + size-1, [](double val) { return val == 0.0; });
    last_one = A[size-1] == 1.0;
    return all_zero && last_one;
}

// Function to extract boundary rows and columns
void ExtractBoundary(double* u0, double* u1, double* u0_top, double* u1_top, double* u0_bottom, double* u1_bottom,
                     double* u0_right, double* u1_right, double* u0_left, double* u1_left,
                     int Nx, int Ny) {
    ExtractRow(u0, u0_top, Nx, Ny - 1);
    ExtractRow(u1, u1_top, Nx, Ny - 1);
    ExtractRow(u0, u0_bottom, Nx, 0);
    ExtractRow(u1, u1_bottom, Nx, 0);
    ExtractColumn(u0, u0_right, Nx, Ny, Nx - 1);
    ExtractColumn(u1, u1_right, Nx, Ny, Nx - 1);
    ExtractColumn(u0, u0_left, Nx, Ny, 0);
    ExtractColumn(u1, u1_left, Nx, Ny, 0);
}

// Perform boundary condition check for top, bottom, right ,and right
void CheckBoundaryConditions(double* u0_top, double* u1_top, double* u0_bottom, double* u1_bottom,
                             double* u0_right, double* u1_right, double* u0_left, double* u1_left,
                             size_t Nx, size_t Ny) {
    BOOST_CHECK_MESSAGE(CheckOneMatrix(u0_top, Nx), "U-velocity Top boundary condition not imposed");
    BOOST_CHECK_MESSAGE(CheckZeroMatrix(u1_top, Nx), "V-velocity Top boundary condition not imposed");
    BOOST_CHECK_MESSAGE(CheckZeroMatrix(u0_bottom, Nx), "U-velocity Bottom boundary condition not imposed");
    BOOST_CHECK_MESSAGE(CheckZeroMatrix(u1_bottom, Nx), "V-velocity Bottom boundary condition not imposed");
    BOOST_CHECK_MESSAGE(CheckZeroOneMatrix(u0_right, Ny), "U-velocity Right boundary condition not imposed");
    BOOST_CHECK_MESSAGE(CheckZeroMatrix(u1_right, Ny), "V-velocity Right boundary condition not imposed");
    BOOST_CHECK_MESSAGE(CheckZeroOneMatrix(u0_left, Ny), "U-velocity Left boundary condition not imposed");
    BOOST_CHECK_MESSAGE(CheckZeroMatrix(u1_left, Ny), "V-velocity Left boundary condition not imposed");
}

/*
// Test case for initial condition (fluid at t=0 is at rest)
BOOST_AUTO_TEST_CASE(BoundaryConditionCheck)
{
    // Domain length for test case
    double Lx_test = 2.5;
    double Ly_test = 1.3;
    // Grid points for test case
    int Nx_test = 19;
    int Ny_test = 7;
    // int Nx_test = 4;
    // int Ny_test = 3;

    double dt_test = 0.003; // time step for test case
    double T_test = 1.2; // Total simulation for test case
    double Re_test = 10; // Reynolds number for test case 
    bool verbose = false; // Do not diplay convergence detail inside this test

    // Create and set up LidDrivenCavity environment
    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetDomainSize(Lx_test,Ly_test);
    solver->SetGridSize(Nx_test,Ny_test);
    solver->SetTimeStep(dt_test);
    solver->SetFinalTime(T_test);
    solver->SetReynoldsNumber(Re_test);
    solver->SetVerbose(verbose);


    // Starting time : t = 0
    solver->Initialise();
    double* u0 = solver->get_u0();
    double* u1 = solver->get_u1();

    // Top and bottom row
    double* u0_top = new double[Nx_test];
    double* u1_top = new double[Nx_test];
    double* u0_bottom = new double[Nx_test];
    double* u1_bottom = new double[Nx_test];
    // Right and left column
    double* u0_right = new double[Ny_test];
    double* u1_right = new double[Ny_test];
    double* u0_left = new double[Ny_test];
    double* u1_left = new double[Ny_test];

    // Extract bottom and top row
    ExtractBoundary(u0, u1, u0_top, u1_top, u0_bottom, u1_bottom, u0_right, u1_right, u0_left, u1_left, Nx_test, Ny_test);
    // Check boundary condition at t=0 (top, bottom, right, and left for u0 and u1) 
    CheckBoundaryConditions(u0_top,u1_top,u0_bottom,u1_bottom,u0_right,u1_right,u0_left,u1_left,Nx_test,Ny_test);


    // Mid point : Integrate half the total simulation time
    double half_time = 0.4;
    solver->IntegrateControl(half_time);
    u0 = solver->get_u0();
    u1 = solver->get_u1();

    // Extract bottom and top row
    ExtractBoundary(u0, u1, u0_top, u1_top, u0_bottom, u1_bottom, u0_right, u1_right, u0_left, u1_left, Nx_test, Ny_test);
    // Check boundary condition at half time (top, bottom, right, and left for u0 and u1) 
    CheckBoundaryConditions(u0_top,u1_top,u0_bottom,u1_bottom,u0_right,u1_right,u0_left,u1_left,Nx_test,Ny_test);


    // End point : Integrate total simulation time
    double end_time = 1.0;
    // Create and set up LidDrivenCavity environment
    LidDrivenCavity* solver2 = new LidDrivenCavity();
    solver2->SetDomainSize(Lx_test,Ly_test);
    solver2->SetGridSize(Nx_test,Ny_test);
    solver2->SetTimeStep(dt_test);
    solver2->SetFinalTime(T_test);
    solver2->SetReynoldsNumber(Re_test);
    solver2->SetVerbose(verbose);

    solver2->Initialise();
    solver2->IntegrateControl(end_time);
    u0 = solver2->get_u0();
    u1 = solver2->get_u1();

    // Extract bottom and top row
    ExtractBoundary(u0, u1, u0_top, u1_top, u0_bottom, u1_bottom, u0_right, u1_right, u0_left, u1_left, Nx_test, Ny_test);
    // Check boundary condition at half time (top, bottom, right, and left for u0 and u1) 
    CheckBoundaryConditions(u0_top,u1_top,u0_bottom,u1_bottom,u0_right,u1_right,u0_left,u1_left,Nx_test,Ny_test);

    BOOST_TEST_MESSAGE("\n  Test results: Boundary condition applied correctly\n");  // Print a message indicating the test case is successful

    delete[] u0;
    delete[] u0_top;
    delete[] u0_bottom;
    delete[] u0_left;
    delete[] u0_right;

    delete[] u1;
    delete[] u1_top;
    delete[] u1_bottom;
    delete[] u1_left;
    delete[] u1_right;
    delete solver2;

}
*/

// Compare two files
bool CompareFiles(const std::string& filename1, const std::string& filename2,int Npts) {
    std::ifstream file1(filename1), file2(filename2);
    std::string line1, line2;

    double tolerance = 1e-10; // absolute error
    double tolerance_percentage = 0.05; // witin 5%
    int tolerance_corrupt = 0;
    bool all_correct = false;
    bool verbose = true;

    // Compare line by line
    while (std::getline(file1, line1) && std::getline(file2, line2)) {
        std::istringstream iss1(line1), iss2(line2);
        double num1, num2;

        // Compare each number in the line
        while (iss1 >> num1 && iss2 >> num2) {

            if (std::abs((num1 - num2)/num2) > tolerance_percentage && std::abs(num1 - num2) > tolerance) {
                if(verbose){
                    std::cout << " Number 1 : "  << num1 << " | Number 2 : " << num2 << " | error : " << std::abs(num1 - num2)<< std::endl;
                }
                // all_correct = false;
                tolerance_corrupt++;
            }
        }

    }

    if(tolerance_corrupt<Npts*0.01){
        all_correct = true;
    }

    return all_correct;
}

// Test case for initial condition (fluid at t=0 is at rest)
BOOST_AUTO_TEST_CASE(FinalResultsCheck)
{
    // Domain length for test case
    double Lx_test = 1;
    double Ly_test = 1;
    // Grid points for test case
    int Nx_test = 9;
    int Ny_test = 9;
    // int Nx_test = 49;
    // int Ny_test = 49;
    int Npts = Nx_test*Ny_test;

    double dt_test = 2e-4; // time step for test case
    double T_test = 1.0; // Total simulation for test case

    double Re_test = 10; // Reynolds number for test case 
    bool verbose = false; // Do not diplay convergence detail inside this test



    // Initialise MPI
    int rank = 0; // ID of process
    int size = 0; // Number of processes
    int err = MPI_Init(NULL,NULL);
    if (err != MPI_SUCCESS) {
        cout << "Error: Failed to initialise MPI" << endl;
    }

    // Get comm rank and size of each process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check if the number of processes is a perfect squares
    int sqrtNum = sqrt(size);
    if(sqrtNum * sqrtNum != size){
		if(rank==0){
			cout << "Please have n^2 of processes" << endl;
		}
		MPI_Finalize();
    }

    // Create a cartesian grid of n x n to compute
    MPI_Comm domain_local;
    int dims[2] = {0, 0};
    int periods[2] = {0, 0};
    int coords[2];
    int reorder = 0; // Do not allow reordering for now

    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &domain_local); // Create the Cartesian communicator
    MPI_Cart_coords(domain_local, rank, 2, coords); // Get the coordinates of the current process in the Cartesian grid

    int rank_left, rank_right, rank_up, rank_down; // Ranks of neighboring processes
    // Get the ranks of the neighboring processes
    // MPI_Cart_shift(domain_local, 0, 1, &rank_up, &rank_down);
    MPI_Cart_shift(domain_local, 0, 1, &rank_down, &rank_up);
    MPI_Cart_shift(domain_local, 1, 1, &rank_left, &rank_right);


    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetDomainSize(Lx_test,Ly_test);
    solver->SetGridSize(Nx_test,Ny_test);
    solver->SetTimeStep(dt_test);
    solver->SetFinalTime(T_test);
    solver->SetReynoldsNumber(Re_test);
    solver->SetVerbose(verbose);
    solver->SetNeighbour(rank_up,rank_down,rank_left,rank_right);
    solver->DomainDecomposition();


    // Lid Driven Cavity Parallel
    std::string output = "test/data/final_test.txt";
    // solver->Initialise();
    // solver->Integrate();


    solver->DomainDecomposition();
    solver->InitialiseParallel();
    solver->IntegrateParallel();
    if(rank==0){
        solver->WriteSolution(output);
        cout << "Successfully written report" << endl;
    }

	MPI_Finalize();

    if(rank==0){
        std::string output_true = "test/data/final.txt";
        // std::string output_true = "test/data/final2.txt";
        bool results = CompareFiles(output,output_true,Npts);
        BOOST_CHECK_MESSAGE(results, "Results do not match"); 
    }

    BOOST_TEST_MESSAGE("\n  Test results: LidDrivenCavity solved correctly\n");

}


BOOST_AUTO_TEST_SUITE_END()

