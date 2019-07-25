These are a sampling of the utility scripts I use in computational 
chemistry workflows. 

casExample.py has an example CASSCF/CASCI/NEVPT2 calculation. The Pyscf package
provides an interface to the Block solver. I use functions and objects from the 
Pyscf module to specify molecule geometry, set the calculation parameters, 
make directories for the calculation to run in, and perform the calculations. 

getEn.py extracts the energies for an arbitrary number of excited states from 
an Orca multiple-structure TD-DFT calculation. The results are printed in a 
table that is easy to parse in other programs. 

LammpsAnalysis/naphDimTraj.py includes a script to calculate three order parameters
for a bridged naphthalene dimer based on a LAMMPS molecular dynamics trajectory. 
Each is plotted in a histogram, which is saved in a pdf file. 

LammpsAnalysis/plotThermo.ipynb is a jupyter notebook demonstrating how to 
make a quick plot of important thermodynamic quantities from a LAMMPS MD run. 

lammpsFileMerge.py merges two LAMMPS format data files (which have info about 
molecule coordinates, bonds, etc). The online LigParGen server can produce 
files for a single molecule, but since I work with dimers, I needed this script
to finish the job. 

makeSandwichDimer.py produces initial coordinates for dimer structures via a 
command-line tool. User can specify center-of-mass offsets for the two molecules 
chosen. By editing the initial coordinate files, the user can flag the script to 
ignore a subset of the atoms. 

pairDist.py calculates the distance between overlapping carbons in a dimer structure. 
Works for any monomer, excludes the covalent linker on the PAH dimer if there is one. 

xyzDimerScan.py sets up several Orca multi-structure TD-DFT calculations to 
examine the intermolecular potential energy surface for a dimer. 
Displacements in the x, y, and z directions are considered. The script produces a 
shell script to submit the jobs to the queue as well. 
