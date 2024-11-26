# 3DCatalystModeling

## Purpose

A general-use 3D solver for coupled gas-phase and aqueous-phase transport/reactions in electroredution catalysts with intermixed hydrophobic and hydrophilic catalysts. The type of catalyst we model/simulate is of the type studied in [this recent paper](https://pubs.acs.org/doi/abs/10.1021/acsenergylett.2c01978). Hydrophobic (porous PTFE) and hydrophilic (porous Cu) domains are interspersed randomly, creating pathways for transport of gaseous and aqueous species, respectively. The Cu phase acts as a catalyst for reduction of CO and CO2 to form simple hydrocarbons, while the PTFE phase helps overcome gas transport limitations in pure Cu catalysts. In this work, we perform geometric restructuring to retain the essential physicochemical dynamicsand scales present in the problem while simplifying the simulation considerably. We are able to recast the problem as a simple unit-problem using a 3D, structured, rectilinear mesh. This page offers a snapshot into the code we use to perform these simulations.

## What is included

This is a brief snapshot of parts of the code, pulled from the active research branch (a private repository belonging to the Mani Group). The intent is to provide a representative sample of my personal approach to scientific software development, as part of a professional portfolio. Here is a brief description of what you will find in each of the main folders/files:

- *source/*: main program files, with a wide variety of different sub-functions (including different versions that were used for experimentation with different schemes)
  - The most important file here is *main.m*, which drives the overall simulation
  - Several separately-defined functions are required to run *main.m*, many of which may be found in this folder
- *MMS_files/*: files for implementation of the *Method of Manufactured Solutions*, a debugging framework in which carefully-designed forcing terms are used to guarantee an analytically known solution, permitting verifiation of the code
- *notes/*: Notes used for documentation and tracking project progress
- *input.json*: Input file, in JSON format. See more details below.
- *launchParallelArray.sh*: The primary way I run simulations is in batch mode on clusters. This is a slurm script lets you launch as many jobs as you would like, to sweep 1D or multi-D paramter spaces.

## What is not included

Since publications related to this code have not been submitted yet, I have elected to withhold all of the post-processing and plotting components. In truth, post-processing and plotting code accounts for as much of the code in the active research repository as the actual simulation components. Additionally, several versions of the *main.m* file were used to prototype and develop different numerical methods and iteration schemes. I have only included the file for the end result of this iterative software design process.

## Running the simulation

I typically run the simulation in batch mode on HPC clusters, permitting quick sweeps of 1D or multi-D paramter spaces. The file to look at is *launchParallelArray.sh*, a SLURM script that can be used to lanch job arrays for sweeping multiple parameters simultaneously. Check out [SLURM documentation](https://slurm.schedmd.com/documentation.html) to learn more about how this file works.

Also be sure to check out *source/runParallelJob.m*, where you can see exactly how parallel jobs are launched. This is also where you can learn details of how to modify *launchParallelArray.sh* in order to sweep multiple paramters simultaneously, instead of just a single parameter (currently, it's a sweep of steady state current densities).

## Input file formatting

The input for the simulation is provided in .JSON format. Here is a quick overview of the different sections of the input file. This will also help you get a quick sense of how flexible and adaptable the software is. **For details about units for each of the parameters, please see the notes/generalNotes.txt file.**

### Directory, time step, operating mode, mesh, and gas parameters

![Input file sample image 1a!](/inputFileImages/input_1a.png "Sample input file, figure 1a.") ![Input file sample image 1b!](/inputFileImages/input_1b.png "Sample input file, figure 1b.")

- Here, you can select the directory where the input file is located and where ouput files will be generated. You can also specify whether the simulation will be restarted from previous output files.
- Select the target time step, in addition to the time-intervals for output generation and plotting (mostly useful for real-time debugging purposes when running on an interactive node). You may choose to sweep the voltage to produce a voltammogram, but I usually run in galvanostatic mode.
- You may ramp the Faradaic reaction kinetic parameters to their final values, if desired: this helps ease the initial transient, improving stability during early times.
- The charge threshhold sets the tolerance for coupling between the equations (poor convergence produces spurious charge, which is non-physical)
- We include the ability to perform equilibration of homogeneous chemistry before the simulation starts, which can be useful for multi-component electrolytes where the initial condition is not in perfect equilibrium. Without this, you might have as rapid transient that causes numerical stability issues during early times.
- Set the electrode potential, and the target current/system capacitance for galavanostatic operation. The external (parallel) capacitor controls the time-scale over which the system approaches the target current.
- This is also where you specify the rectilinear domain size for the aqueous region, and the mesh size in each dimension.
  - Mesh refinement near the boundaries is an option here. The "LRC" variable lets you specify if you would like refinement (down to the "min" parameter) near the left boundary, right boundary, or both boundaries.
- Finally, this is also where all of the gas-phase parameters are set, including composition, domain size, and flow rate. Porosity and tortuosity of the gas-transporting region can be set appropriately, in addition to pressure of the gas supply.

### Boundary conditions and physical constants

![Input file sample image 2!](/inputFileImages/input_2.png "Sample input file, figure 2.")

- Boundary conditions in each of the dimensions (no-flux at hydrophilic/hydrophobic interfaces)
- The "y Right" boundary is the interface between the catalyst and bulk electrolyte
  - Dirichlet conditions are used, corresponding to the known electrolyte properties
  - The electric potential is fixed
- Physical consants are assigned standard values

### Porous structure and aqueous species

![Input file sample image 3!](/inputFileImages/input_3.png "Sample input file, figure 3.")

- There is only a single "layer" here (the "layer" terminology is a vestige of antoher project for simulation of multi-layered cells, which uses a similar input file format)
- Properties of the porous structure may be set here
- The aqueous species, diffusion coefficients, valences, and initial concentrations are set here
- You can specify the LaTeX names, for fancy plot labeling

### Homogeneous reactions

![Input file sample image 4!](/inputFileImages/input_4.png "Sample input file, figure 4.")\

- Elementary homogeneous reactions and rates
- The Wien coefficient is a vestigal parameter from a previous study on bipolar membranes, it is typically not used here
  - But, the option exists, in case it should become relevant in a future application

### Faradaic reactions

![Input file sample image 5!](/inputFileImages/input_5.png "Sample input file, figure 5.")

- Electron-transferring reaction at the pore walls in the Cu phase
- **For details about units for each of the parameters, please see the notes/generalNotes.txt file.**

## Mesh generation

Users may use exponential mesh refinement near one or both boundaries in each direction. Specification of a mimimum and maximum mesh size is enough to analytically solve for the correct analytical stretching transformation that yields a cell distribution with smooth exponential transition from the smallest to largest cells. Users may select refinement near a single boundary in particular or both boundaries in each of the dimensions.

## Boundary conditions

Users have the ability to set no-flux or Dirichlet boundary conditions in each of the dimensions. Physically, we expect a no-flux condition for aqueous species at interfaces between aqueous and gaseous regions. At the interface between the catalyst and the outer bulk electrolyte, we set Dirichlet conditions for the species concentrations and the electric potential, corresponding to the known electrolyte composition and the applied voltage.

## Numerical methods

We offer a brief summary of the key numerical methods; an exhaustive description will appear in an upcoming publication. Second-order central differences are used for spatial discretization, and a third-order, implicit, multi-step method is used for temporal discretization. The start-up steps use lower-order implicit schemes in time, since the third-order method requires data from three previous points in time (which are not available for clean starts). Crucially, we use an iterative approach to coupling the equations. Furthermore, the selection of implicit numerics in the narrowest dimension improves numerical stability while enabling parallel decomposition of the domain. Thus, examining *main.m*, you will find that the **par-for** loop is used to parallelize the aqueous-phase solve for each time step.

## Output format and processing

Separate binary output files are created for the aqueous region and gaseous region, labeled with the simulation name and the time step. I have provided functions for reading in the binary files and producing human-readable MatLab data structures, which can be found in *source/*

## Notes files

- **Please see the notes/generalNotes.txt file for details about units and definitions of input parameters**
- The *notes/doNext.txt* file is a track record of objectives set and completed in the active research branch
