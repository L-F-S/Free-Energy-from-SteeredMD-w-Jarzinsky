# Nonequilibrium Thermodynamic Integration by means of Steered Molecular Dynamics of Equus Caballus Acylophosphatase (1APS)

This repository contains scripts to calculate a free energy profile of [1APS](https://www.rcsb.org/structure/1APS "1APS pdb") between its native and denatured state, by applying Jarzynski's identity \cite{jarzynski1997nonequilibrium} to a dataset of steered molecular dynamics simulations, starting from a PDB file of the coordinates of the native state.
The theory behind the project is described in more detail [here](https://it.overleaf.com/read/ygvjwnpbbrfk "Nonequilibrium Thermodynamic Integration by means of Steered Molecular Dynamics of Equus Caballus Acylophosphatase (1APS)").
![](img/aps_unfolded.jpg)
![](img/aps_folded.jpeg)


# Directories
```bash
├── steeredMDsimulations
│   ├── 
│   ├── input_files/
│   │  ├── aps.pdb
│   │  ├── aps.
├── Jarzinsky equality calculation
│   ├── 
├── README.md
```
# Requirements:
  * [VMD](www.ks.uiuc.edu/Research/vmd), a molecular graphics program
  * [NAMD](www.ks.uiuc.edu/Research/namd), a molecular dynamics simulation program
  * Python 3
# 1. Run SMD simulations

See [namd tutorial guide on steered molecular dynamics simulations](https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node16.html) to know how to set up and run constant velocity pulling steered MD simulations. You should end up with somAll files are already created inside  /steeredMDsimulaitons/input_files/  The configuration file for NAMD is in /steeredMDsimulaitons/apf.conf

## 1.1 Input:
-[protein pdb](https://www.rcsb.org/structure/1APS "1APS pdb")
-protein.psf file  
## Command:

## Output:

# 2. Analyze data and calculate free energy:


# 3. Plot what you like
![](img/W_wrt_dist.png)
