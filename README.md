 The lid-driven cavity problem is a classical benchmark in computational fluid dynamics that 
 provides valuable insights into the behaviour and efficiency of the CFD algorithm. 
 In this problem, a square cavity is filled with fluid, and the motion of the fluid is driven 
 by the continuous movement of the top lid, while the other boundaries remain stationary.
 

To run main solver:
1.  make
    ./solver

or

(You can specify the solver arguments in make file without using explicity typing ./solver)
2. make run Np=4 Nx=801 Ny=801 T=5


To run test:
1.  make unittests
    ./testing (Please wait for a few minutes if using 1 process)

or

2. make run_test Np=9 Nt=2


To generate HTML and LATEX from doxygen (please configure file from the makefile guidelines (echo) if generate Doxyfile)
    cd docs/
    doxygen


To visualise simulation results:
Open the FlowVisualisation.ipynb (change results file path if necessary)
Run all and view results from the output



