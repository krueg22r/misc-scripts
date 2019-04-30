#! /usr/local/bin/python3

import argparse

# takes coords from file, puts in conventient array format
# also returns number of atoms 
def procXYZ(fileName): 
	xyzFile = open(fileName, 'r') 
	lines = xyzFile.readlines()
	nAt = int(lines[0])
	coordArray = list()
	for line in lines[2:]: 
		tokens = line.split()
		coordArray.append(tokens)
	return nAt, coordArray

# get center of mass for given coord array, excluding: 
# hydrogens and carbons that have exclude symbol as token 

def getCom(coords): 
	massDict = {'C':12.00, 'O':16.0, 'N': 14.0}
	com = [0, 0, 0]
	totMass = 0
	for line in coords: 
		if len(line) == 4: #exclude ones w/extra flag
			elem = line[0]
			if elem != 'H': 
				mass = 0
				try: 
					mass = massDict[elem]
				except IndexError: 
					print('element '+elem+' not found! ') 
					exit()
				for i in range(0,3): 
					com[i] += (float(line[i+1]) * mass)
				totMass += mass 
	for i in range(0,3): 
		com[i] = com[i]/totMass
	return com 
		

def main(): 
	parser = argparse.ArgumentParser()
	parser.add_argument("xyz1", type=str, help = "first xyz file")
	parser.add_argument("xyz2", type=str, help = "second xyz file")
	parser.add_argument("-x", type = float, default = 0.0)
	parser.add_argument("-y", type = float, default = 0.0)
	parser.add_argument("-z", type = float, default = 3.3)
	args = parser.parse_args()
	comOffset = [float(args.x), float(args.y), float(args.z)]
	nAt1, bigCoords = procXYZ(args.xyz1) 
	nAt2, smallCoords = procXYZ(args.xyz2)
	if nAt1 < nAt2: 
		temp = bigCoords
		bigCoords = smallCoords 
		smallCoords = temp 
	bigCom = getCom(bigCoords) 
	smallCom = getCom(smallCoords) 
	deltaCom = [0,0,0]
	for i in range(0,3): 
		deltaCom[i] = smallCom[i] - bigCom[i]
	print(str(nAt1+nAt2)+'\n')
	for line in smallCoords: 
		print(line[0], end = '   ') 
		for i in range(1,4): 	
			#print('{:>13.10}'.format(str(float(line[i])-deltaCom[i-1]+comOffset[i-1])), end = '')
			print('{0:0.7f}'.format(float(line[i])-deltaCom[i-1]+comOffset[i-1]), end = '  ')
		print('')
	for line in bigCoords: 
		print(line[0], end = '   ')
		for i in range(1,4):
			print('{:>13.10}'.format(line[i]), end = '')
		print('')
	
main()
