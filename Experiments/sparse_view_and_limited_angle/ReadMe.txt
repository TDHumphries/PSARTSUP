ReadMe: Code for Superiorized pSART
Author: T. Humphries
Last updated: Nov 9 2017

Experiments directory

This directory contains scripts for the numerical experiments presented in "Superiorized algorithm for reconstruction of CT images from sparse-view and limited-angle polyenergetic data" by Humphries, Winn and Faridani (https://arxiv.org/abs/1701.03396) The scripts are:

i)   make_FORBILD_data.m: Generates data from the FORBILD head phantom which is used in all experiments. So, run this one first. It will save a file called FORBILD_data.mat in the directory.
ii)   make_FORBILD_data_80kVp.m: Same as above but generates data with a 80 kVp spectrum, used as an additional set of experiments in the paper.
iii)  exp1_sparseview.m: 			 Runs sparse-view experiments.
iv) exp1_limitedangle.m: 			 Runs limited-angle experiments.

