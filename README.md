# USP (User-controlled Scatter Plots) for SPINS
Extended mixing diagnostics in MATLAB for SPINS

## Introduction to method (what is USP, why do we want to use it?)

USP (User-controlled Scatter Plots), or joint probability plots tell us where two variables overlap, so we can find out if, when and where we get combinations of two variables. 

USP is a statistical method (it's Baysean P(A|B)) of analysing a physical domain by decomposing into a 2-dimensional paired histogram of two fluid properties. In this case, one property is usually a physical property of the fluid (density, temperature, salinity, tracer concentration), and the other a measure of flow (speed, Kinetic Energy, vorticity, enstrophy, dissipation). 

Typically, the sorting algorithm can be used to investigate mixing, telling how much of each type of energy there is at each timestep, and how energy moves between those types (see [plot\_diagnos.m](https://github.com/ddeepwel/SPINSmatlab/blob/master/plotting/plot_diagnos.m) in SPINSmatlab). Unlike other methods typically used to investigate mixing, it is primarily qualitative, but it informs us where mixing takes place, and the fluid parcels that have undergone mixing, as well as the fact we can change the physical domain of interest easily (which we can't do for a sorting algorithm). 

![Schematic of USP](./F1_SchematicQSP.png)

## Research questions available using this method
- Tracing how fluid parcels change as they move through physical space
- Tagging of fluid parcels - we can map back from USP to physical space
- Understanding covariance of fluid properties

## Features \& Usage
<details>
<summary>
  USP_2D/USP_3D
</summary>

[2D MATLAB Tool](usp_2d.m)

[3D MATLAB Tool](usp_3d.m)

MATLAB tool for calculating the USP histograms for 2D//3D simulations, mapped or unmapped grids (easy to adjust parameters and run local analysis) - mapped only for 2D. 

Both USP\_2d, USP\_3d have the same input variables (wich broadly match their equivalents in USP_to_physical):
- ii is the output number,
- var1 is the first variable (USP x axis),
- var2 second variable (USP y axis),
- spat\_lims is the spatial region to plot \[xmin xmax zmin zmax\]
- var\_lims is the limits of the variables (essentially axis limits of USP), NOTE this is \[var2min var2max var1min var2max\], and can be set only for var2. It's probably usually good practice to set the var\_lims, as SPINS often outputs one or two super off the scale values which just drown out any other data (this is just like we usually would set the caxis).
- doPlot is a boolean switch to make plots, or to output data. 

MATLAB files are set up to use as variable 1, density, salinity, or any other SPINS direct output file. As variable 2, it can do KE, speed, enstrophy, vorticity, dissipation, or any other SPINS direct output file. 

</details>

<details>
<summary>
  USP_to_physical
</summary>
  
[2D MATLAB Tool](usp_to_physical.m)

[3D MATLAB Tool](usp_to_physical_3d.m)

MATLAB tool to identify what fluid in physical space meets the criteria identified from a USP histogram.

In interactive mode (Region variable not input) USP\_to\_physical (2D only) will present to you the output from USP\_mapped (it in fact simply calls that function), prompting you to select the Region of Interest in the USP histogram (by clicking two opposite corners of the rectangle). It will then plot the fluid within this region of interest for both variables and the USP OR output the boolean mask for the ROI.

- ii is the output number,
- var1 is the first variable (USP x axis),
- var2 second variable (USP y axis),
- spat\_lims is the spatial region to plot \[xmin xmax zmin zmax\]
- var\_lims is the limits of the variables (essentially axis limits of USP), NOTE this is \[var2min var2max var1min var2max\], and can be set only for var2. It's probably usually good practice to set the var\_lims, as SPINS often outputs one or two super off the scale values which just drown out any other data (this is just like we usually would set the caxis).
- Region is the Region of Interest to be plotted

For the following command: 
`>> usp_to_physical(67, 'rho', 'KE', [5.5 7], [0 0.005 -0.0095 0.0095]);`

A region of interest can be isolated interactively:
![](./usp_to_physical.gif)


MATLAB files are set up to use as variable 1, density, salinity, or any other SPINS direct output file. As variable 2, it can do KE, speed, enstrophy, vorticity, dissipation, or any other SPINS direct output file. N.B.  when we use a square measure you do essentially squash low values together, and stretch the high values togehter, so consider this. 
There's also an "isInvert" switch in qsp\_to\_physical which needs to be manually changed in the code, if it's set to true it shows us the physical space which is outside the highlighted usp space - in some ways it's "What part of the flow isn't interesting"

TODO: Add a schematic here
TODO: qsp\_to\_physical to work with the csv read in from SPINSqsp

</details>

<details>
<summary>
    SPINS tools
</summary>
  
- [SPINS qsp](https://git.uwaterloo.ca/SPINS/SPINS_main/-/tree/master/src/cases/qsp) - A c++ tool for calculating USP within the SPINS architecture. Best if using large outputs, and/or 3D simulations. 
- [UNDER DEVELOPMENT] SPINS qsp read tool - Tool to read the .csv output by the SPINS qsp
- [UNDER DEVELOPMENT] MATLAB QSP Statistics tool. 
</details>


### Choosing the criterion for USP space
1. Based on the initial conditions (focus on the region that has changed)
2. An interesting feature on the USP space
3. Some kind of pre-determined criteria (exceeding a threshold)
Perhaps may be a bit hand-wavey (we know that super high KE values are probably bs, low KE is boring). 

## See also:
- [SPINS\_main](https://git.uwaterloo.ca/SPINS/SPINS_main.git)
- [SPINSmatlab](https://github.com/ddeepwell/SPINSmatlab.git)


## References
Penney et al., 2020; Diapycnal mixing of passive tracers by Kelvin-Helmholtz instabilities. JFM, [https://doi.org/10.1017/jfm.2020.483](https://doi.org/10.1017/jfm.2020.483)
