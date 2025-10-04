# FallingAxisymmetricObject

[![Falling objects](figure/FallingObjectSimulation.png)](https://youtu.be/4_bd3tomnOg)

This repository contains the simulation and 3D printing files used in the *Falling Object* demonstration for the course **MT2461-25** at TU Delft, lectured by Prof. Gabriel Weymouth.  
- Physics demo: Manuel Cabral  
- Simulation scripts: Tzu-Yao Huang  

Slides for the class are available [here](https://manuel-cabral.github.io/Hydromechanica_tutorial/).


## Repository Structure

- **ThreeD_Object/** - STL files for 3D printing (objects are sliced in half for printing).  
- **scripts/** - Julia and Python scripts for running and post-processing simulations, and theoretical calculations.  
- **figure/** - Example results, plots, and demo images.


## Running the Simulation (Julia)

1. **Install Julia**  
   Follow the [official instructions](https://julialang.org/install/) or run:
   ```bash
   curl -fsSL https://install.julialang.org | sh
   ```

2. **Clone this repository**

   ```bash
   git clone https://github.com/TzuYaoHuang/FallingAxisymmetricObject.git
   cd FallingAxisymmetricObject/scripts
   ```

3. **Set up the environment**

   ```bash
   julia --project -E 'using Pkg; Pkg.add(url="https://github.com/weymouth/BiotSavartBCs.jl.git"); Pkg.instantiate()'
   ```

4. **Run the simulation**

   ```bash
   julia --project -t auto Falling.jl all 96
   ```

   * Replace `all` with:

     * `run` → only run the simulation
     * `pp` → only post-process (plot the drag coefficient)
   * `96` is the grid resolution in lateral directions (can be changed but should be multiplier of 16).

5. **View results**

   * Output `.pvd` files can be opened in **ParaView**.
   * Drag coefficient comparison plots and the data are in the `figure/` folder.


### Numerical methods and packages

This simulation demos uses a few Julia packages to make the computational domain manageable, runnable on GPU, and to flexibly define shapes:

* **WaterLily.jl**: a lightweight, high-performance incompressible flow solver in Julia, which supports multi-threaded CPU or GPU backends and immersed boundaries. 
* **BiotSavartBCs.jl**: you can shrink the computational domain (i.e. use a smaller grid) while still resolving key flow features, since the boundary influence (vortex contributions) is handled more “analytically” at domain edges through Bio-Savart kernel.
* **ParametricBodies.jl**: a companion package to WaterLily for defining smooth parametric bodies (2D / 3D) via parametric curve representations, where NURBS is also possible.
  * We use ParametricBodies.jl to define the teardrop shape used in the falling object simulation. It provides `ParametricBody` types (a subtype of `WaterLily.AbstractBody`) that define the signed distance, normals, and motion, enabling custom shapes without manually writing complex geometry functions.
