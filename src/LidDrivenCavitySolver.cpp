#include <iostream>
#include <stdio.h>
#include <chrono>
#include <ctime>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "LidDrivenCavity.h"


/**
 * @brief Entry point for solving the 2D lid-driven cavity incompressible flow problem.
 * 
 * The main file sets up MPI, creates a Cartesian grid,
 * initializes the solver, integrates the solution, and writes the final solution to a file.
 * It also measures and displays the time taken for the computation.
 * 
 * @param argc Number of command-line arguments.
 * @param argv Program arguments to specify simulation configuration
 * @return An integer exit status
 */
int main(int argc, char **argv)
{
	// Program arguments
	double Lx; // Length of the domain in the x-direction
	double Ly; // Length of the domain in the y-direction
	double dt; // Time step 
	double T; // Total simulation time
	double Re; // Reynolds number in computation
	int Nx; // Number of grid points in x-direction
	int Ny; // Number of grid points in y-direction
	bool verbose = true; // Display more hint message
	bool legacy = false; // Execute legacy baseline numerical code
    int root = 0; // Root process
    std::string folder_results = "results/";



    // Extracting program arguments with BOOST library
    po::options_description opts(
        "Solver for the 2D lid-driven cavity incompressible flow problem");
    opts.add_options()
        ("Lx", 	po::value<double>(&Lx)->default_value(1.0), "Length of the domain in the x-direction")
        ("Ly", 	po::value<double>(&Ly)->default_value(1.0), "Length of the domain in the y-direction")
        ("dt", 	po::value<double>(&dt)->default_value(2e-4), "Time step ")
        ("T", 	po::value<double>(&T)->default_value(1.0), "Total simulation time")
        ("Re", 	po::value<double>(&Re)->default_value(10), "Reynolds number used in computation")
        ("Nx",	po::value<int>(&Nx)->default_value(9), "Number of grid points in x-direction")
        ("Ny", 	po::value<int>(&Ny)->default_value(9), "Number of grid points in y-direction")
        ("verbose", po::value<bool>(&verbose)->default_value(true), "Display more hint message")
        ("legacy", po::value<bool>(&legacy)->default_value(false), "Execute legacy serial baseline numerical solver")
		("help", "Display argument details");
	// Ensure arguments are parsed correctly
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, opts), vm);
        po::notify(vm);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
	
	// Display argument details and terminates program if help flag
	if(vm.count("help")){
		std::cout << opts << std::endl;
		return 0;
	}
	


    // Initialise MPI
    int rank = 0; // ID of process
    int size = 0; // Number of processes
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        cout << "Error: Failed to initialise MPI" << endl;
        return -1;
    }

    // Get comm rank and size of each process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check if the number of processes is a perfect squares
    int sqrtNum = sqrt(size);
    if(sqrtNum * sqrtNum != size){
		if(rank==root){
			cout << "Please have n^2 of processes" << endl;
		}
		MPI_Finalize();
		return 0;
    }

    // Check domain length
    if (Lx <= 0 || Ly <= 0 || Nx <= 0 || Ny <= 0) {
        if (rank == root) {
            std::cout << "Please ensure positive domain length and valid number of nodes." << std::endl;
        }
        MPI_Finalize();
        return 0;
    }

    // Allow only positive time arguments
    if(T<=0 || dt<=0){
		if(rank==root){
			cout << "Please have a positive time argument for the solver" << endl;
		} 
		MPI_Finalize();
		return 0;    
    }

    // Allow only positive Reynolds number
    if(Re<=0){
		if(rank==root){
			cout << "Please have a positive Reynolds number for the solver" << endl;
		} 
		MPI_Finalize();
		return 0;    
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
    MPI_Cart_shift(domain_local, 0, 1, &rank_down, &rank_up);
    MPI_Cart_shift(domain_local, 1, 1, &rank_left, &rank_right);


    // Retrieving the number of threads
    int Nt =1;
    #pragma omp parallel
    {
        Nt = omp_get_num_threads();
    }
    // Prevent the usage of too many threads (especially when unspecified) which slows down coarse mesh problems
    if (Nt > 25) {
        if(rank==0){
            std::cout << "Warning: Too many threads (" << Nt << ") to be used, limiting the number of threads to 1\n" << std::endl;
        }
        omp_set_num_threads(1);
        Nt = 1;
    }


    // Begin program and displays current time
    auto start = std::chrono::system_clock::now();
    if(rank==root){
        std::time_t time_now = std::chrono::system_clock::to_time_t(start);
        std::cout << "2D lid-driven cavity incompressible flow problem | Time :  " << std::ctime(&time_now);
    }

    // Core program
    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetDomainSize(vm["Lx"].as<double>(), vm["Ly"].as<double>());
    solver->SetGridSize(vm["Nx"].as<int>(),vm["Ny"].as<int>());
    solver->SetTimeStep(vm["dt"].as<double>());
    solver->SetFinalTime(vm["T"].as<double>());
    solver->SetReynoldsNumber(vm["Re"].as<double>());
    solver->SetVerbose(verbose);
    solver->SetThreads(Nt);
    solver->SetNeighbour(rank_up,rank_down,rank_left,rank_right);
    solver->DomainDecomposition();

    // Print fluid solver configuration if root rank
    if(rank==root){
        solver->PrintConfiguration();
    }

    // Warn about legacy mode if enabled and root rank
    if (legacy && rank==root) {
        std::cout << "Warning: Legacy mode is enabled. Executing legacy baseline numerical code.\n" << std::endl;
    }

    // Initialize solver
    if (not(legacy)) {
        solver->InitialiseParallel();
    } else {
        solver->Initialise();
    }
    
    // Write initial solution if root rank
    if(rank==root){
        solver->WriteSolution(folder_results+"ic.txt");
    }

    // Integrate solver
    if (not(legacy)) {
        solver->IntegrateParallel();
    } else {
        solver->Integrate();
    }

    // Write final solution if root rank
    if(rank==root){
        solver->WriteSolution(folder_results+"final.txt");
    }

	// End of Program ad displays current time
    if(rank==root){
        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "Time spent: " << elapsed_seconds.count() << " seconds" << std::endl;
        std::cout << "End of program | Time : " << std::ctime(&end_time);
    }



    // Finalise MPI.
    MPI_Finalize();
	return 0;
}
