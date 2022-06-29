# QSP for SPINS

## Introduction to method (what is QSP, why do we want to use it?)

QSP, or joint probability plots tell us where two variables overlap, so we can find out if, when and where we get combinations of two variables. 

QSP is a statistical method of analysing a physical domain by decomposing into a 2-dimensional paired histogram of two fluid properties. In this case, one property is usually a physical property of the fluid (density, temperature, salinity, tracer concentration), and the other a measure of flow (speed, Kinetic Energy, vorticity, enstrophy, dissipation). 

It is a Baysean technique that finds the unity of the two properties. 

Unlike other methods typically used to investigate mixing, it is primarily qualitative, but it informs us where mixing takes place, and the fluid parcels that have undergone mixing, as well as the fact we can change the physical domain of interest easily (which we can't do for a sorting algorithm). 


## Research questions available using this method
- Tracing how fluid parcels change as they move through physical space
- Tagging of fluid parcels - we can map back from QSP to physical space

## Tools provided:
- MATLAB QSP_mapped (easy to adjust parameters and run local analysis)
- MATLAB QSP_to_physical
- SPINS qsp (best for large simulations or 3D)
- SPINS qsp read tool
- MATLAB Statistics tool. 

MATLAB files are set up to use as variable 1, density, salinity, or any other SPINS direct output file. As variable 2, it can do KE, speed, enstrophy, vorticity, dissipation. 

There's a switch in qsp_to_physical (isInvert) which needs to be manually changed in the code, if it's set to true it shows us the physical space which is outside the highlighted qsp space - in some ways it's "What part of the flow isn't interesting"

Note, when we use a square measure you do essentially squash low values together, and stretch the high values togehter, so consider this. 

## Choosing the criterion for QSP space
1. Based on the initial conditions (focus on the region that has changed)
2. An interesting feature on the QSP space
3. Some kind of pre-determined criteria (exceeding a threshold)
Perhaps may be a bit hand-wavey (we know that super high KE values are probably bs, low KE is boring). 

## References
Penney et al 2020
