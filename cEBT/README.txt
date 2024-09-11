README for cEBT

This folder contains files related to the implementation of our connectivity generalization to the Entropic Bonding Theory

Descriptions of files are:

1). sample_binary_cEBT.m
		Main script for computing cEBT orbital and bonding orbitals between patchy, connected NPs

2). generate_meshgrid.m
		Convenient function for meshgrid generation for cEBT calculations

3). generate_kernel_evaluator.m
		Convenient function for parameterization the surface of an anisotropic particle for cEBT calculations

4). compute_wavefunction_core.m
		Computes the shape orbitals from cEBT for the core (shaped) NPs

5). compute_wavefunction_bond.m
		Compute the shape orbitals from cEBT for the patches on NPs

6). compute_energy.m
		Compute the total cEBT energy from the computed wavefunction for each NP

7). calc_rscale.m
		Convenient function to compute the location of highest pP density from the computed wavefunctions 	