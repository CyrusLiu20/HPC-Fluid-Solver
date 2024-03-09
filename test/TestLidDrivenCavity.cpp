#include "../src/LidDrivenCavity.h"
#include <iostream>

#define BOOST_TEST_MODULE LidDrivenCavityClass
#include <boost/test/included/unit_test.hpp>

// Test case for initial boundary condition value
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
    double Dx = solver->get_dx();
    double Dy = solver->get_dy();
    BOOST_CHECK_EQUAL(Dx,Lx_test/(Nx_test-1));
    BOOST_CHECK_EQUAL(Dy,Ly_test/(Ny_test-1));
}

// // Test case for initial boundary condition value
// BOOST_AUTO_TEST_CASE(InitialBoundaryValueCheck)
// {

// }