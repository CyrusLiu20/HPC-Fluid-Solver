#include <iostream>
#include <boost/program_options.hpp>
using namespace std;

namespace po = boost::program_options;

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
	bool verbose; // Display more hint message
	std::string program_description = "Solver for the 2D lid-driven cavity incompressible flow problem:";
	
    // Define options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("Lx", 	po::value<double>(&Lx)->default_value(1.0), "Length of the domain in the x-direction")
        ("Ly", 	po::value<double>(&Ly)->default_value(1.0), "Length of the domain in the y-direction")
        ("dt", 	po::value<double>(&dt)->default_value(2e-4), "Time step ")
        ("T", 	po::value<double>(&T)->default_value(1.0), "Total simulation time")
        ("Re", 	po::value<double>(&Re)->default_value(10), "Reynolds number used in computation")
        ("Nx",	po::value<int>(&Nx)->default_value(9), "Number of grid points in x-direction")
        ("Ny", 	po::value<int>(&Ny)->default_value(9), "Number of grid points in y-direction")
        ("verbose", po::value<bool>(&verbose)->default_value(true), "Display more hint message")
		("help", "Display argument details");


	// Ensure arguments are parsed correctly
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
	
	
	// Display argument details and terminates program if help flag
	if(vm.count("help")){
		std::cout << program_description << std::endl;
		std::cout << desc << std::endl;
		return 0;
	}
	
		
	return 0;
}


