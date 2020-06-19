# WRF_LES_diagnostics
modification of the WRF model source code (version 4.1) to derive turbulent quantities and mean flow diagnostics.

Replace files in the directories dyn_em and Registry with the files inlcuded in this repository. 

Diagnostics include terms extracted from the dynamical core of the WRF model and explicitly derived terms.

Resolved variances of wind components can be calculated in post-processing by:

- var(U) = UUMEAN - UMEAN * UMEAN

- var(V) = VVMEAN - VMEAN * VMEAN

- var(W) = WWMEAN - WMEAN * WMEAN

Diagnostics include average potential temperature tendencies (K s^-1) resulting from:

- subgrid-scale (SGS )turbulent mixing in each cartesian coordiante (SGSTRBx, SGSTRBy, SGSTRBz)

- microphysics and radiation parameterization (T_TENDM_MP and T_TENDM_RAD)

- advection with the mean flow diagnosed from first-order differences in cartesian coordiantes (MADVx, MADVy, MADVz)

- resolved part of turbulent mixing (RESTRBx, RESTRBy, RESTRBz)

Cartesian components of the resolved part of turbulent mixing can be calculated in postprocessing by:

- RESTRBx = meanTOTADVx - MADVx

- RESTRBy = meanTOTADVy - MADVy

- RESTRBz = meanTOTADVz - MADVz

Averageing is based on the default WRF code in module_avgflx.F. Activate averaging by setting do_avgflx_em = 1 in the namelist section &dynamics. Averaging period is equal to the primary output interval of the simulation (namelist variable history_interval). Output of individual quantities can be activated either by directly modifying registry.avgflx or by using the "Run-time IO" capability of the WRF model (https://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap5.htm#runtimeio)
