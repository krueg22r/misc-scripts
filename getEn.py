#! /usr/bin/env python3

# A program to extract excited-state energies from an Orca multiple structure output file
# It is used for dimers in which molecule separations vary in only one dimension. 
# Works for output from any kind of TDDFT calculation (even thought the outputs are not the same) 

import argparse 

def main():
	parser = argparse.ArgumentParser()
	# this is the Orca output file with an excited-state multiple structure scan
	parser.add_argument('outfile', type=str)
	args = parser.parse_args()
	outName = args.outfile
	out = open(outName, 'r')
	eGround = 0
	deltaEs = list() # excitation energies 
	outfileLines = out.readlines()
	rCoord = 1
	for i, line in enumerate(outfileLines): 
		if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
			line1 = i
			tokens = outfileLines[line1+2].split()
			rVal = float(tokens[rCoord])
			print('{:>8.4}'.format(rVal), end='')	
		if ' RI-MP2 CORRELATION ENERGY' in line:
			tokens = line.split()
			eGround += float(tokens[3])
		if 'Total Energy       :' in line: 
			tokens = line.split()
			eGround += float(tokens[3])
		if 'STATE  ' in line and 'au' in line: 
			tokens = line.split()
			deltaEs.append(float(tokens[3]))
		if 'Dispersion correction' in line: 
			tokens = line.split()
			eGround += float(tokens[2])
		if 'FINAL SINGLE POINT ENERGY' in line: 
			# print final energies 
			for en in deltaEs: 
				 print('{:>14.9}'.format(eGround+en), end='')
			print('')
			deltaEs = list()
			eGround = 0
			
				
main()
