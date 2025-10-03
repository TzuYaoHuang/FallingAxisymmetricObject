# FallingAxisymmetricObject

[![Falling objects](figure/FallingObjectSimulation.png)](https://youtu.be/4_bd3tomnOg)

This repository contains the simulation and 3D printing files used in the *Falling Object* demonstration for the course **MT2461-25** at TU Delft, lectured by Prof. Gabriel Weymouth.  
- Physics demo: Manuel Cabral  
- Simulation scripts: Tzu-Yao Huang  

Slides for the class are available [here](https://manuel-cabral.github.io/Hydromechanica_tutorial/).

---

## Repository Structure

- **ThreeD_Object/** - STL files for 3D printing (objects are sliced in half for printing).  
- **scripts/** - Julia and Python scripts for running and post-processing simulations, and theoretical calculations.  
- **figure/** - Example results, plots, and demo images.

---

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
   julia --project -E 'using Pkg; Pkg.add(url="https://github.com/weymouth/BiotSavartBCs.jl.git", rev=true); Pkg.instantiate()'
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

