/*
Three Test cases for LidDrivenCavity class

Program arguments check
1. DomainCheck: Check correct mesh size
2. ReynoldsCheck: Check correct Reynolds number and Nusselt number specified

Numerical simulation check
3. InitialConditionCheck: Check if initial condition (v,s,u0,u1) is applied correctly
*/

#include "../src/LidDrivenCavity.h"
#include <iostream>
#include <boost/algorithm/cxx11/all_of.hpp>

#define BOOST_TEST_MODULE LidDrivenCavityClass
#include <boost/test/included/unit_test.hpp>

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

    delete solver;
}