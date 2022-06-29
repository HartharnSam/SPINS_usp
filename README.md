# QSP for SPINS

## Introduction to method (what is QSP, why do we want to use it?)

QSP, or joint probability plots tell us where two variables overlap, so we can find out if, when and where we get combinations of two variables. 

QSP is a statistical method (it's Baysean P(A|B) ) of analysing a physical domain by decomposing into a 2-dimensional paired histogram of two fluid properties. In this case, one property is usually a physical property of the fluid (density, temperature, salinity, tracer concentration), and the other a measure of flow (speed, Kinetic Energy, vorticity, enstrophy, dissipation). 

Typically, the sorting algorithm can be used to investigate mixing, telling how much of each type of energy there is at each timestep, and how energy moves between those types (see [plot\_diagnos.m](https://github.com/ddeepwel/SPINSmatlab/blob/master/plotting/plot_diagnos.m) in SPINSmatlab). Unlike other methods typically used to investigate mixing, it is primarily qualitative, but it informs us where mixing takes place, and the fluid parcels that have undergone mixing, as well as the fact we can change the physical domain of interest easily (which we can't do for a sorting algorithm). 

## Research questions available using this method
- Tracing how fluid parcels change as they move through physical space
- Tagging of fluid parcels - we can map back from QSP to physical space

## Tools provided:
- [QSP\_mapped](qsp_mapped.m) - MATLAB tool for calculating the QSP histograms (easy to adjust parameters and run local analysis)
- [QSP\_to\_physical](qsp_to_physical.m) - MATLAB tool to identify what fluid in physical space meets the criteria identified from a QSP histogram. 
- \*SPINS qsp (best for large simulations or 3D)
- \*SPINS qsp read tool
- \*MATLAB Statistics tool. 
\* under development/not included here

MATLAB files are set up to use as variable 1, density, salinity, or any other SPINS direct output file. As variable 2, it can do KE, speed, enstrophy, vorticity, dissipation, or any other SPINS direct output file. 

### Usage of tools
Both QSP\_mapped and QSP\_to\_physical have the same input variables, ii is the output number, var1 is the first variable, var2 second variable, xlims is the x region to plot, and ke_lims is the limits of the second variable. %TODO: Flag to sort out both xlims and ke\_lims to work on two both axes - 2x2 matrix kind of thing.

QSP\to\_physical will then present to you the output from QSP\_mapped (it in fact simply calls that function), prompting you to select the space in the QSP histogram (by clicking two opposite corners of the rectangle). The interactive functionality can be switched off and replaced by fixed values by switching around some commented lines. 

There's also an "isInvert" switch in qsp\_to\_physical which needs to be manually changed in the code, if it's set to true it shows us the physical space which is outside the highlighted qsp space - in some ways it's "What part of the flow isn't interesting"

It's probably usually good practice to set the ke\_lims, as SPINS often outputs one or two super off the scale values which just drown out any other data (this is just like we usually would set the caxis). 

N.B. Note, when we use a square measure you do essentially squash low values together, and stretch the high values togehter, so consider this. 

### Choosing the criterion for QSP space
1. Based on the initial conditions (focus on the region that has changed)
2. An interesting feature on the QSP space
3. Some kind of pre-determined criteria (exceeding a threshold)
Perhaps may be a bit hand-wavey (we know that super high KE values are probably bs, low KE is boring). 

## References
Penney et al., 2020; Diapycnal mixing of passive tracers by Kelvin-Helmholtz instabilities. JFM, [https://doi.org/10.1017/jfm.2020.483](https://doi.org/10.1017/jfm.2020.483)
