ReadMe: Code for Superiorized pSART
Author: T. Humphries
Last updated: Nov 9, 2017

Reconstruction directory

Contains the following files:

i) pSART_sup.m: main method for iterative reconstruction. Allows user to perform monoenergetic reconstruction (SART), polyenergetic (pSART), and choose to superiorize the method with a secondary objective of their choice. Also allows for weighted least squares option as described in Humphries and Gibali. See experiments folder for examples of how to call the method.

ii) get_*_params.m, print_params.m: helper methods used to read input to pSART_sup.m and print the options back out to the user.

iii) create_sysmats.m: generates system matrix and other matrices needed for iterative reconstruction using the MIRT. Automatically opens up a parallel pool of workers to help accelerate this process.

iv)  project_poly.m: implements polyenergetic forward projection

v) grad_TV.m, grad_ATV.m, compute_omega.m: used to compute the TV and ATV functions for superiorization.