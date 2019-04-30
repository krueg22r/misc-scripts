#! /usr/bin/env python

from pyscf import gto
from pyscf import scf
from pyscf.mcscf import avas
from pyscf import tools
from pyscf import mcscf
from pyscf import mrpt
from pyscf import dmrgscf
import time
import os

# Intermolecular separation for this run 
x = 3.1
# Set up directories for run 
rootDir = '/panfs/panbox/home/rkrueger/Heterodimer/ButadieneBenz/Pyscf/cc-pvtz/S1/'+str(x)+'/'
runName = rootDir+'Runtime'
scratchName = rootDir+'Scratch'
outputName = 'out.pyscf_'+str(x)
os.system('mkdir '+runName)
os.system('mkdir '+scratchName)
os.system('date > startTime')
# Settings for DMRG solver 
dmrgscf.settings.BLOCKSCRATCHDIR=scratchName
dmrgscf.settings.BLOCKRUNTIMEDIR=runName
dmrgscf.settings.MPIPREFIX ='/panfs/panbox/home/rkrueger/intel-2013/impi/4.1.3.048/intel64/bin/mpirun -n 8'
# Basis set 
basisName = 'cc-pvtz'
# Build molecule
mol = gto.Mole()
mol.build(
	verbose=4,
atom = [['C', (0.682701, 1.753076, x)], ['C', (0.534938, 0.428158, x)], ['C', (-0.726008, -0.321109, x)], ['C', (-1.960308, 0.182633, x)], ['H', (1.665845, 2.204406, x)], ['H', (-0.162196, 2.431298, x)], ['H', (1.431809, -0.183803, x)], ['H', (-0.617379, -1.401412, x)], ['H', (-2.152068, 1.248981, x)], ['H', (-2.826760, -0.465085, x)], ['C', (-0.705820, 2.463432, 0.0)], ['C', (0.498589, 1.768065, 0.0)], ['C', (0.498589, 0.377333, 0.0)], ['C', (-0.705820, -0.318033, 0.0)], ['C', (-1.910229, 0.377333, 0.0)], ['C', (-1.910229, 1.768065, 0.0)], ['H', (-0.705820, 3.545494, 0.0)], ['H', (1.435682, 2.309097, 0.0)], ['H', (1.435682, -0.163698, 0.0)], ['H', (-0.705820, -1.400095, 0.0)], ['H', (-2.847323, -0.163698, 0.0)], ['H', (-2.847323, 2.309097, 0.0)]],
	basis=basisName,
	output=outputName,
	max_memory=16000
	)
# Hartree-Fock calculation 
mf = scf.RHF(mol)
mf.kernel()
mf.analyze()
# Get active space for CASSCF calculation
n_orb_active_space, nelec_act_space, new_all_orbs = avas.kernel(mf, 'C 2pz')
# CASSCF 
mc = dmrgscf.DMRGSCF(mf, n_orb_active_space, nelec_act_space)
mc.state_specific_(state=1)
mc.fcisolver.maxM = 500
mc.kernel(new_all_orbs)
mc_orbs = mc.mo_coeff
# CASCI 
mc = mcscf.CASCI(mf, n_orb_active_space, nelec_act_space) 
mc.fcisolver = dmrgscf.DMRGCI(mol)
mc.fcisolver.maxM = 1200
mc.fcisolver.nroots = 2
mc.kernel(mc_orbs)
# NEVPT2 correction 
mps_nevpt_e2 = mrpt.NEVPT(mc,root=1).compress_approx(maxM=1200).kernel()
mc.analyze()
os.system('date > endTime')
