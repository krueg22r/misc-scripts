#! /panfs/panbox/home/rkrueger/python/bin/python3

# A program to compute several order parameters for the bridged
# naphthalene dimer based on a LAMMPS MD trajectory. 
# Parameters distributions are plotted with histograms. 

import argparse
import numpy as np
import matplotlib.pyplot as plt
import math

def getOrderParams(coords, nAtoms): 
	# indices of "overlapping" carbons in the eclipsed structure 
        ccList = [[3,41],[2,39],[6,37],[4,36],[10,35],[8,33],[12,31],[11,29],[16,28],[14,25]] 
        nPairs = 10
        meanParam  = 0
        maxParam = 0
	# calculate distance between overlapping carbons 
        for pair in ccList: 
                dist = np.linalg.norm(coords[pair[0]-1]-coords[pair[1]-1])
                meanParam += dist
                maxParam = max(dist, maxParam) 
        meanParam /= nPairs
	# calculate angles between molecule and linker carbons
        vec1a = coords[16]-coords[21]
        vec2a = coords[40]-coords[21]
        vec1b = coords[21]-coords[16]
        vec2b = coords[2]-coords[16]
        dot1 = np.dot(vec1a,vec2a)
        dot2 = np.dot(vec1b,vec2b)
        theta1 = math.acos(dot1 / (np.linalg.norm(vec1a)*np.linalg.norm(vec2a)))*180/math.pi
        theta2 = math.acos(dot2 / (np.linalg.norm(vec1b)*np.linalg.norm(vec2b)))*180/math.pi
        return meanParam, maxParam, theta1, theta2

# make histograms of order parameters 
def makePlot(datArray, xLabel, yLabel, fileName):
        plt.figure(figsize=(6,6))
        fontSize = 12
        plt.rcParams.update({'font.size': fontSize})
        plt.xlabel(xLabel)
        plt.ylabel(yLabel) 
        n, bins, patches = plt.hist(datArray, bins=50, normed = True, facecolor='b')
        plt.savefig(fileName+'.pdf', format='pdf', bbox_inches='tight')

def main():
        parser = argparse.ArgumentParser()
        parser.add_argument('outfile', type=str)
        args = parser.parse_args()
        outName = args.outfile
        out = open(outName, 'r')

        line = out.readline()
        frameCount = 0
        nAtoms = 0
        timestep = 0.000001*500 # units of ns 
        skipFrames = 1/timestep
        maxParams = list() # order parameter based on maximum distance between paired carbons
        meanParams = list() # order parameter based on mean distance between paired carbons
        angles = list()
        minSampleVal = 5
        maxSampleVal = 11
        nBins = 6
        delta = (maxSampleVal-minSampleVal)/nBins
        sampleCount = np.zeros((nBins))
        while True:
                if line == '': 
                        break
                if frameCount < skipFrames: 
                        continue
                if 'NUMBER OF ATOMS' in line: 
                        line = out.readline()
                        nAtoms = int(line)
                        frameCount += 1
                if 'ATOMS id x y z' in line: 
			# read coordinates 
                        coords = np.zeros((nAtoms, 3))
                        for i in range(0, nAtoms): 
                                line = out.readline()
                                tokens = line.split()
                                for j in range(1,4): 
                                        coords[i,j-1] = float(tokens[j])
                        meanParam, maxParam, angle1, angle2 = getOrderParams(coords, nAtoms)
                        maxParams.append(maxParam)
                        meanParams.append(meanParam) 
                        angles.append(angle1)
                        angles.append(angle2)
                        print('meanParam = '+str(meanParam))
                        bin = (meanParam - minSampleVal)//delta
                        print('bin = '+str(bin))
                line = out.readline()
        meanParams = np.asarray(meanParams)
        maxParams = np.asarray(maxParams)
        angles = np.asarray(angles) 
	# Plot resutls 
        makePlot(meanParams, 'Mean C-C distance ($\mathrm{\AA}$)', 'Probability', 'meanParam') 
        makePlot(maxParams, 'Max C-C distance ($\mathrm{\AA}$)', 'Probability', 'maxParam') 
        makePlot(angles, 'Linker-PAH angle ($^{\circ}$)', 'Probability', 'angle') 
        

main()
