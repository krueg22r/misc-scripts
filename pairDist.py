#! /usr/bin/env python3

# finds the distance between nearest-neighbor carbons in 
# two different PAH groups. Excludes the linker if one is present. 

import glob
import math
import numpy as np
from itertools import combinations 
from inertia import arrayInertia

# run in directory with output files 
outfiles = glob.glob('out.complex*')
for file in outfiles:
        out = open(file, 'r') 
        outLines = out.readlines()
        # determine whether job finished
        if 'TERMINATED' not in outLines[len(outLines)-2]:
                continue #if not, skip this file
        nAtom = 0
        lineNum = 0
        calcFinished = False 
        for index, line in enumerate(outLines):
                if 'Number of atoms' in line:
                        tokens = line.split()
                        nAtom = int(tokens[4])
        # everything printed after this is about the final, optimized struct
                if('FINAL ENERGY' in line):
                        calcFinished = True
                        lineNum = index+6
                # calculate mass, read in coords of heavy atoms 
                        cCoords = list()
                        hCoords = list()
                        oxContaining = '0'
                        for j in range(lineNum,lineNum+nAtom):
                                tokens = outLines[j].split()
                                elem = tokens[0]
                                if elem == 'H': 
                                        xyz = [float(i) for i in tokens[1:]]
                                        hCoords.append(xyz)
                                else: 
                                        xyz = [float(i) for i in tokens[1:]]
                                        cCoords.append(xyz) 
                        nCs = len(cCoords) 
                        nH = np.zeros(nCs)
			# Figure out how many hydrogen atoms each carbon has 
                        for h in hCoords: 
                                neigh = 0; 
                                minDist = 100; 
                                for j in range(0, nCs): 
                                        hCoords = np.array(h) 
                                        cCoordArr = np.array(cCoords[j])
                                        dist =  np.linalg.norm(hCoords-cCoordArr)
                                        if dist < minDist: 
                                                minDist = dist
                                                neigh = j
                                nH[neigh] += 1 
                        molCoords = list()
                        linkerCoords = list()
			# if a carbon has more than one H, it's in the linker 
                        for j in range(0, nCs):
                                if nH[j] > 1: 
                                        linkerCoords.append(cCoords[j])
                                else: 
                                        molCoords.append(cCoords[j])
                        # also assign atoms to molecules 
                        mol1indices = set() 
                        mol2indices = set() 
                        pairSet = set()
                        mol1indices.add(0) # place the first atom here 
                        coords = np.array(molCoords) # we have now eliminated linkers Cs! 
                        nHeavyAt = len(coords) 
                        # Determine which atoms are in which monomer 
                        bondCounter = 0; 
                        for i in range(0, nHeavyAt-1): 
                                for j in range(i+1, nHeavyAt): 
                                        # monomer membership 
                                        bondCutoff = 1.7
                                        dist = np.linalg.norm(coords[i]-coords[j])
                                        if dist < bondCutoff: 
                                                pairSet.add((i, j)) 
                                                bondCounter += 1; 
                        done = False 
			# Add first set of connected edges to mol 1 
                        while (not done): 
                                pairAdded = (-1, -1) 
                                for pair in pairSet: 
                                        added = False
                                        if (pair[0] in mol1indices and pair[1] not in mol1indices):
                                                mol1indices.add(pair[1])
                                                added = True
                                        if (pair[1] in mol1indices and pair[0] not in mol1indices):
                                                mol1indices.add(pair[0])
                                                added = True
                                        if added: 
                                                pairAdded = pair 
                                                break 
                                if pairAdded[0] ==  -1: #we failed to add one pair tho we tried them all 
                                        done = True
                                else: 
                                        pairSet.remove(pairAdded) 
                        # the leftovers go to mol2 
                        for pair in pairSet: 
                                if pair[0] not in mol2indices and pair[0] not in mol1indices: 
                                        mol2indices.add(pair[0])
                                if pair[1] not in mol2indices and pair[1] not in mol1indices: 
                                        mol2indices.add(pair[1])
                        bigMonIndices = mol1indices
                        smallMonIndices = mol2indices
                        if len(mol2indices) > len(mol1indices): 
                                bigMonIndices = mol2indices
                                smallMonIndices        = mol1indices        
                        coordsBig = list()
                        coordsSmall = list()
                        for ind in bigMonIndices: 
                                coordsBig.append(coords[ind])
                        for ind in smallMonIndices: 
                                coordsSmall.append(coords[ind])
                        coordsBig = np.array(coordsBig)
                        coordsSmall = np.array(coordsSmall) 
                        bigCom = np.average(coordsBig, axis=0) 
                        smallCom = np.average(coordsSmall, axis=0)
                        com = np.average(np.concatenate((coordsBig,coordsSmall),axis=0),axis=0)
                        bigComDists = list()
                        smallComDists = list()
                        ccDists = list()
                        bigLinkerC = -1
                        smallLikerC = -1
                        smallMinDist = 100
                        bigMinDist = 100
                        counter = 0
                        for atom in linkerCoords: 
                                distToBig = np.linalg.norm(bigCom-atom)
                                distToSmall = np.linalg.norm(smallCom - atom) 
                                if distToBig < bigMinDist: 
                                        bigMinDist = distToBig
                                        bigLinkerC = counter
                                if distToSmall < smallMinDist: 
                                        smallMinDist = distToSmall
                                        smallLinkerC = counter
                                counter += 1
                        for i in range(0, len(coordsBig)): 
                                for j in range(0, len(coordsSmall)): 
					# molecules mostly lie in z-plane 
                                        xyThresh = 1.0
                                        xyDist = np.linalg.norm(coordsBig[i,:2]-coordsSmall[j,:2])
                                        if xyDist < xyThresh: # overlapping criterion 
                                                dist = np.linalg.norm(coordsBig[i]-coordsSmall[j])
                                                comDistBig = np.linalg.norm(coordsBig[i]-com)
                                                comDistSmall = np.linalg.norm(coordsSmall[j]-com)
                                                bigComDists.append(comDistBig)
                                                smallComDists.append(comDistSmall) 
                                                ccDists.append(dist)
                        # sort in order of dist from center of mass, print
                        smallComDists, ccDists = zip(*sorted(zip(smallComDists, ccDists)))
			# print resutls 
                        print('# atom from center of mass,      intermolecular distance (\AA)')
                        for i in range(0, len(ccDists)): 
                                print(str(i+1)+'  ', end='')
                                print(format(ccDists[i], '.4f'))
                       	# done with this output file  
                        break
