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

## Running the simulation

## Input file formatting

## Mesh generation

## Numerical methods

## Method of Manufactured Solutions

## Notes files

##
