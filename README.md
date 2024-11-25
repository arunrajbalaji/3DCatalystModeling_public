# 3DCatalystModeling

## Purpose

A general-use 3D solver for coupled gas-phase and aqueous-phase transport/reactions in electroredution catalysts with intermixed hydrophobic and hydrophilic catalysts. The type of catalyst we model/simulate is of the type studied in [this recent paper](https://pubs.acs.org/doi/abs/10.1021/acsenergylett.2c01978). Hydrophobic (porous PTFE) and hydrophilic (porous Cu) domains are interspersed randomly, creating pathways for transport of gaseous and aqueous species, respectively. The Cu phase acts as a catalyst for reduction of CO and CO2 to form simple hydrocarbons, while the PTFE phase helps overcome gas transport limitations in pure Cu catalysts. In this work, we perform geometric restructuring to retain the essential physicochemical dynamicsand scales present in the problem while simplifying the simulation considerably. We are able to recast the problem as a simple unit-problem using a 3D, structured, rectilinear mesh. This page offers a snapshot into the code we use to perform these simulations.

## What is included

This is a brief snapshot of parts of the code, pulled from the active research branch (a private repository belonging to the Mani Group). The intent is to provide a representative sample of my personal approach to scientific software development, as part of a professional portfolio. Here is a brief description of what you will find in each of the main folders/files:

- *source/*: main program files, with a wide variety of different sub-functions (including different versions that were used for experimentation with different schemes)
  - The most important file here is *main.m*, which drives the overall simulation
  - Several separately-defined functions are required to run *main.m*, all of which may be found in this folder
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

The input for the simulation is provided in .JSON format. Here is a quick overview of the different sections of the input file. This will also help you get a quick sense of how flexible and adaptable the software is.

### Directory, time step, and mesh parameters

Here, you can select the directory where the input file is located and where ouput files will be generated. You can also specify whether the simulation will be restarted from previous output files.

Select the target time step, in addition to the time-intervals for output generation and plotting (mostly useful for real-time debugging purposes when running on an interactive node). You may choose to sweep the voltage to produce a voltammogram, but I usually run in galvanostatic mode. You may ramp the Faradaic reaction kinetic parameters to their final values, if desired: this helps ease the initial transient, improving stability during early times.

We include the ability to perform equilibration of homogeneous chemistry before the simulation starts, which can be useful for multi-component electrolytes where the initial condition is not in perfect equilibrium. Without this, you might have as rapid transient that causes numerical stability issues during early times.

Set the electrode potential, and the target current/system capacitance for galavanostatic operation. The external (parallel) capacitor basically controls the time-scale over which the system approaches the target current.

This is also where you specify the rectilinear domain size for the aqueous region, and the mesh size in each dimension. Mesh refinement near the boundaries is an option here. The "LRC" variable lets you specify if you would like refinement (down to the "min" parameter) near the left boundary, right boundary, or both boundaries.

Finally, this is also where all of the gas-phase parameters are set, including composition, domain size, and flow rate. Porosity and tortuosity of the gas-transporting region can be set appropriately, in addition to pressure of the gas supply.

![Input file sample image 1!](/inputFileImages/input_1.png "Sample input file, figure 1.")

Boundary condition selection occurs here. 

![Input file sample image 2!](/inputFileImages/input_2.png "Sample input file, figure 2.")

## Mesh generation

## Boundary conditions

## Numerical methods

## Output format and processing

## Method of Manufactured Solutions

## Notes files

##
