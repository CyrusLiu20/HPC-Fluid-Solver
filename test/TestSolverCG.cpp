/*
Test cases for LidDrivenCavity class

Solver Check
1. Arbitrary matrix check: Vorticity and stream function as arbitrary matrix
2. Real case check: real and realistic v and s taken from LidDrivenCavity class
*/

#include "../src/SolverCG.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/test/unit_test.hpp>

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


BOOST_AUTO_TEST_SUITE_END()

