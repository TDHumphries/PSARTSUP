# PSARTSUP
Code for reconstruction of polyenergetic CT data using superiorization

This archive contains code used for the Superiorized pSART method described in "Superiorized algorithm for reconstruction of CT images from sparse-view and limited-angle polyenergetic data" by Humphries, Winn and Faridani (https://arxiv.org/abs/1701.03396), and "Superiorized polyenergetic reconstruction algorithm for reduction of metal artifacts in CT images" by Humphries & Gibali (http://faculty.washington.edu/thumphri/papers/HG17.pdf). I am releasing it under the GNU Public License: https://www.gnu.org/licenses/gpl-3.0.en.html

The "core" directory contains functions implementing the reconstruction algorithms, routines for generating projection and phantom data, and implementations of the TV and anisotropic TV calculations. The "Experiments" directory contains scripts for running the experiments described in the aforementioned papers. Both directories contain ReadMes describing their contents.

The code uses the Matlab Image Processing toolbox, as well as the Michigan Image Reconstruction Toolbox (MIRT) to generate polyenergetic X-ray spectra and attenuation coefficients for different materials, as well as construct system matrices. You can find the MIRT at https://web.eecs.umich.edu/~fessler/code/ . Note that you may need to modify the file irt/ct/xray_read_dens.m included in that toolbox as it does not contain density information for some of the body tissues used in the experiments. (These can be found in the online NIST database at http://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html). 

You are free to use this code for educational purposes, for your own experiments or to compare against other methods. I would ask that you cite the following paper in any publications which use this code:

i) T. Humphries, J. Winn and A. Faridani. Superiorized algorithm for reconstruction of CT images from sparse-view and limited-angle polyenergetic data. Phys. Med. Biol. 62(16), p. 6762. 2017. http://iopscience.iop.org/article/10.1088/1361-6560/aa7c2d

Additionally, if using for experiments involving metal objects or using the weighted least squares option, please cite:

ii) T. Humphries and A. Gibali. "Superiorized polyenergetic reconstruction algorithm for reduction of metal artifacts in CT images." 2017 IEEE Nuclear Science Symposium Conference Record.

The code also draws on work presented in the following papers:

"An efficient polyenergetic SART (pSART) reconstruction algorithm for quantitative myocardial CT perfusion" by Lin & Samei (2014) http://dx.doi.org/10.1118/1.4863481
"Superiorization: An optimization heuristic for medical physics" by Herman, Garduno, Davidi and Censor (2012) http://dx.doi.org/10.1118/1.4745566
"A limited-angle CT reconstruction method based on anisotropic TV minimization" by Chen, Jin, Li and Wang (2013) https://doi.org/10.1088/0031-9155/58/7/2119
"Simulation tools for two-dimensional experiments in x-ray computed tomography using the FORBILD head phantom" by Yu, Noo, Dennerlein, Wunderlich, Lauritsch and Hornegger. https://doi.org/10.1088/0031-9155/57/13/N237

TH
