Readme for MD Simulations 
Paper: "Elucidating the Interplay Between Entropy-Driven and Patch-Mediated 
Bonding in Directing Nanoscale Assemblies"
Authors: Kireeti Akkunuri, Xiangyu Zhang, Thi Vo
Dated: September 2024

------------------------------------------
Part 1
Filename: 'Part1_Main_Initialization.m'
This contains MATLAB code to initialize shape chains, which outputs configurational text files.
Associated function files: 'rot_quat.m','rot_az.m','gen_poly.m','calc_inertial_tensor.m'
To set the shape: Assign the relevant shape ('Tetrahedron', 'Octahedron', or 'Cube') to the variable 'shape_name'. 
The associated input shape templates are included as TXT files.
To set the chain length: Assign the value to 'N'. 
(The values of N used for the study are: 40, 60, 80, 100, 120, 140, 160, 180, 200)
To set the patch-bonding location: assign 'corner', 'edge', or 'face' to bond_location.
Output: 'init.txt', 'init_bond.txt', 'verts.txt', 'verts_rigid.txt'

--------------------------------------------
Part 2
Filename: Part2_MD_run.py
This contains the run script for HOOMD-Blue (v2) in python. 
It takes as input the configurational files generated from Part 1. 
The documentation on how to use HOOMD-Blue can be found at: https://hoomd-blue.readthedocs.io/en/v2.9.7/
To run multiple simulations for statistics, the random number 'seed' variable has to be changed.
Output: GSD trajectory file
(https://gsd.readthedocs.io/en/v2.6.1/)

--------------------------------------------
Part 3
Filename: Part3_Analysis.py
This reads the GSD trajectory to calculate the size. 
It utilizes the Freud library to simplify box unwrapping calculations. 
(https://freud.readthedocs.io/en/v2.12.1/)
A sample GSD file ('sample_single_part1.gsd') is attached for convenience. 
This corresponds to a corner-bonded tetrahedron of chain length 40. 
The output is in the form of time-averaged radius of gyration, which is then averaged across 3 different runs.
Output: 'r_g.txt' (Not attached) 
