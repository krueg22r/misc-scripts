#! /usr/bin/env python3

# Sets up a potential energy surface scan for the benzene-naphthalene 
# dimer. Coordinates are provided for differently-rotated benzene molecules. 
# This script scans x, y, and z coordinates. 

import os
import glob
import numpy as np
from scriptFunc import writeScriptLong

def printCoords(naphLines, benzCoords, elem, outFile, dx, dy, dz): 
	dr = [dx, dy, dz]
	print("30", file=outFile)
	print("", file=outFile)
	for i in range (0,12): 
		print(elem[i], end = '  ', file=outFile)
		for q in range(0,3): 
			print('{:>16.8}'.format(benzCoords[i,q]+dr[q]), end='', file=outFile)
		print('', file=outFile)
	for i in range(2, len(naphLines)): 
		print(naphLines[i], end = '', file=outFile)


def main():
	benzFiles = glob.glob('benz*')
	naphFiles = glob.glob('naph*')
	naphFileName = naphFiles[0]
	naphFile = open(naphFileName, 'r')
	naphLines = naphFile.readlines()
	naphCoords = np.zeros((6,3))
	for i in range(6, 12): #consider right-hand C ring
		tokens = naphLines[i].split()
		for q in range(0,3):
			naphCoords[i-6,q] = tokens[q+1]
	naphCom = np.mean(naphCoords, axis=0)
	#now loop over angular displacements
	angleCounter = -1
	for b in benzFiles: 
		angleCounter += 1
		#first tranform benz coords so benz sits in in the middle
		benzFile = open(b, 'r') 
		benzLines = benzFile.readlines()
		benzCoords = np.zeros((12,3))
		benzElems = list()
		for i in range(2,14): 
			tokens = benzLines[i].split()
			benzElems.append(tokens[0])
			for q in range(0,3): 
				benzCoords[i-2,q] = tokens[q+1]
		benzCom = np.mean(benzCoords, axis=0) 
		for i in range(0,12): 
			for q in range(0,3): 
				benzCoords[i,q] += (naphCom[q]-benzCom[q])
				if q == 1: 
					benzCoords[i,q] += 1.2124
				if q == 2: 
					benzCoords[i,q] += 2.8
		xMin = 0
		xMax = 1.4
		yMin = 0
		yMax = 2.4
		zMin = 2.8
		zMax = 3.3
		baseName = "angle"+str(angleCounter)		
		dirName = baseName+"scan"
		os.system("mkdir "+dirName)
		xyzFile = open(baseName+".xyz", 'w')
		vmdFile = open("vmd"+baseName+str(angleCounter)+".xyz", 'w')
		orcaFile = open(baseName, 'w') 
		#print job info to orca file
		print("! RKS RI-JK D3BJ B2PLYP def2-TZVP def2/JK TightSCF Grid5 NoFinalGrid", file=orcaFile) 
		print("", file=orcaFile)
		print("%MaxCore 1800", file=orcaFile)
		print("%pal nprocs 12", file=orcaFile)
		print("end", file=orcaFile)
		print("", file=orcaFile)
		print("%tddft", file=orcaFile)
		print("NRoots 1", file=orcaFile)
		print("MaxDim 8", file=orcaFile)
		print("end", file=orcaFile)
		print(" ",  file=orcaFile)
		print("* xyzfile 0 1 "+baseName+".xyz", file=orcaFile)
		print("", file=orcaFile)
		writeScriptLong(baseName, 12)
		os.system("mv "+baseName+" "+dirName)
		os.system("mv "+baseName+".sh "+dirName)
		dx = 0.2
		dy = 0.2
		dz = 0.1
		nx = int((xMax-xMin)/dx) +2
		ny = int((yMax-yMin)/dy) + 2
		nz = int((zMax-zMin)/dz) + 1
		for z in range(0,nz): 
			for y in range(0,ny): 
				for x in range(0,nx): 
					printCoords(naphLines, benzCoords, benzElems, vmdFile, x*dx, y*dy, z*dz)
					printCoords(naphLines, benzCoords, benzElems, xyzFile, x*dx, y*dy, z*dz)
					print(str(round(x*dx,2))+" "+str(round(y*dy,2))+" "+str(round(z*dz,2)))
					#if z != (nz-1) and x != (nx-1) and y != (ny-1): 
					print(">", file=xyzFile)
		#print("", file=xyzFile)
		os.system("mv "+"vmd"+baseName+str(angleCounter)+".xyz"+" "+dirName)
		os.system("mv "+baseName+".xyz "+dirName)	
		
		
	
	
	
main()
