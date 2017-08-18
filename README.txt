This package contains the code and input data needed to reproduce the analyses of Larsen, Ma, and Kaplan, "How important is connectivity for surface-water fluxes? A generalized expression for flux through heterogeneous landscapes." Explanations of all code and inputs follows.

This package requires Matlab, along with its Statistics Toolbox and Optimization Toolbox, to run.

A. INPUT FILES
1. Synthetic Landscapes input only.xlsx: This file contains the results of the processed SWIFT2D simulations, with no column headings, for ease of reading into Matlab. Metadata associated with this file is as follows: 

Column A: Water depth (m) at which simulation was run
Column B: Effective water depth (m) for landscapes with patch coverage = 0.1
Column C: Effective water depth (m) for landscapes with patch coverage = 0.35
Column D: Effective water depth (m) for landscapes with patch coverage = 0.425
Column E: Effective water depth (m) for landscapes with patch coverage = 0.5
Column F: Effective water depth (m) for landscapes with patch coverage = 0.575
Column G: Effective water depth (m) for landscapes with patch coverage = 0.65
Column H: Effective water depth (m) for landscapes with patch coverage = 0.9

Note that "effective water depth" is the same as Heff in the companion paper. The file is 840 rows long. The first 210 rows correspond to 30 landscape realizations with an anisotropy of 6. Rows 211-420 correspond to 30 landscape realizations with an anisotropy of 4. Rows 421-630 correspond to 30 landscape realizations with an anisotropy of 2. Rows 631-840 correspond to 30 landscape realizations with an anisotropy of 1. 

2. Synthetic Landscapes Q.xls: Contains the "raw" results of the SWIFT2D simulations. These are discharge values, in m^3/s. Metadata is as follows: 

Column A: Discharge (m^3/s) for landscapes with patch coverage = 0.1
Column B: Discharge (m^3/s) for landscapes with patch coverage = 0.35
Column C: Discharge (m^3/s) for landscapes with patch coverage = 0.425
Column D: Discharge (m^3/s) for landscapes with patch coverage = 0.5
Column E: Discharge (m^3/s) for landscapes with patch coverage = 0.575
Column F: Discharge (m^3/s) for landscapes with patch coverage = 0.65
Column G: Discharge (m^3/s) for landscapes with patch coverage = 0.9

3. Landscape metrics.xls: This file contains all of the landscape metric data used in the model for omega. Metadata is as follows:
Column A: Name of the synthetic landscape. The number that follows the leading "R" is the patch coverage, as a percent. Each patch coverage has 120 domains. The first 30 have an anisotropy of 6, the second 30 have an anisotropy of 4, the third 30 have an anisotropy of 2, and the final 30 have an anisotropy of 1.
Column B: Weighted, directed DCI for patches
Column C: Weighted, undirected DCI for patches
Column D: Weighted, directed DCI for channels
Column E: Weighted, undirected DCI for channels
Column F: Unweighted,  directed DCI for channels
Column G: Unweighted, undirected DCI for channels
Column H: fd: Patch fractal dimension. Varies between 1 and 2, with 2 being more Euclidean (i.e., less complex edges)
Column I: Mean spanning path tortuosity, directed
Column J: Mean spanning path tortuosity, undirected


B. ANALYSIS CODES
1. OmegaModelingLowFlow.m: This code reads in simulation data, calculates omega, and explores functional forms for a low-flow model of omega. The first part of the code explores the functional form of omega at h = zp. The second part explores the functional form of omega at h<zp.
2. OmegaModelingHighFlow.m: This code reads in simulation data, calculates omega,and explores functional forms for a high-flow model of omega, starting with the slope and ending with the intercept. The final portion of the code checks the overall omega fit and then the velocity fit.
3. CalculateFlowSensitivity.m: This code conducts a generalized sensitivity analysis given a model of omega. The final part of the code addresses the case-study: it determines how much of a change in velocity was due to change in each landscape configuration factor over a trajectory of landscape degradation. 
4. CalculateLowFlowSensitivity.m: This code conducts a generalized sensitivity analysis given a model of omega. In the beginning of the code, the user selects between two options: h = zp, or h < zp. 
5. GenFig3B.m: This routine performs a comparison between a field-calibrated Kadlec model of flow through Everglades site U3 (by Choi and Harvey, Wetlands, 2014; hereafter abbreviated CH) an upscaled Manning model with optimized parameters, and a standard Manning model (through a homogeneous environment) with optimized roughness. It then generates the plot depicted as Fig. 3B.
6. OmegaTheoretical.m: This routine generates computations of omega for theoretical solutions of flow over a landscape with patches in parallel and a landscape with patches in series. This code was used to generate Figs. 1A-B. 