# WRF_LES_diagnostics
modification of the WRF model source code (version 4.1) to derive turbulence and mean flow diagnostics 

Replace files in the directories dyn_em and Registry with the files inlcuded in this repository. 

Diagnostics include terms extracted from the dynamical core of the WRF model and explicitly derived terms.

Resolved variances of wind components can be calculated in post-processing by:

var(U) = UUMEAN - UMEAN * UMEAN

var(V) = VVMEAN - VMEAN * VMEAN

var(W) = WWMEAN - WMEAN * WMEAN

Cartesian components of the resolved part of turbulent mixing can be calculated in postprocessing by:

RESTRBx = meanTOTADVx - MADVx

RESTRBy = meanTOTADVy - MADVy

RESTRBz = meanTOTADVz - MADVz
