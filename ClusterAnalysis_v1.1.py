#!/usr/bin/python

# To Do list:
# 
# polish the code.

oneline = "Cluster analyze dump file and write to PDB file."

docstr =  """ 
Cluster analysis used to study polyelectrolyte complexation
Read in : dump file generated by LAMMPS

Atom type 1 : Neutral bead in polycation (PC)
Atom type 2 : Charged bead in polycation   (+PC)
Atom type 3 : Neutral bead in polyanion  (PA)
Atom type 4 : Charged bead in polyanion   (-PA)
Atom type 5 : negatively charged counterions for polycation chains (-PC counterion)
Atom type 6 : positively charged counterions for polyanion chains (+PA counterion)
Atom type 7 : negative ions for added salt (- salt)
Atom type 8 : positive ions for added salt (+ salt)

input variable explanation:

counterion_size    : size of counterion, varies from 0.1 to 1.0 here is 0.5
len_pan            : length of polyanion chain, 128
len_pcat           : length of polycation charin, 32
num_pan            : number of polyanion, 32 
charge_ratio       : x = Cc/Ca, where x is charge_ration, 
                     Cc is number of beads in polycation, 
                     Ca is number of beads in polyanion
salt_concentration : s = Cs/Ca, where s is salt_concentration, 
                     Cs is number of salt bead in solution,
                     
Usage: (This is for Shengfeng's run data file)
$ ./ClusterAnalysis.py dumpfile len_pan len_pcat num_pan charge_ratio salt_concentration

Example:
../binary_IPEC/tot_charge_ratio_2.0/la128_lc32_r2_npa32_s2/dump_atom_pe_la128_lc32_r2_npa32_s2_equi.10000000

$ ./ClusterAnalysis.py dump_atom_pe_la128_lc32_r2_npa32_s1_equi_x0.5.20000000  128 32 32 0.5 1.0

charge_fraction = 0.5 is constant in proposal.

Author: Yanfei Tang & Dr. Shengfeng Cheng 02/04/2016

Usage of v1.1:
$ ./ClusterAnalysis.py dumpfile len_pan num_pan len_pcat num_pcat salt_concentration

Example:
$ ./ClusterAnalysis.py dump.pe.la128.na32.lc32.nc80.rho0.01.r2.s0.2.equi.10000000 128 32 32 80 0.2
"""

# -----Library import -----
import dump as dp
import numpy as np
import sys
import time
import string
import os

# -----Read From Terminal Input ------
filename = sys.argv[1]
len_pan = int(sys.argv[2])
num_pan = int(sys.argv[3])
len_pcat = int(sys.argv[4])
num_pcat = int(sys.argv[5])
salt_concentration = float(sys.argv[6])


# class definition

class ClusterAnalysis:
    
    """
    Atributes:
    ----------
    ifile : string
        name of input dump file. 
    
    beadSize : float64
        size of bead in PA/PC polymer chain.
    counterionSize : float64
        size of counterion.
    saltSize : float64
        size of salt atom.
    distFactor : float
        the factor that controls how close two mols needed to be in order to belong to the same cluster.
    numCounterion : int
        total number of counterion
    numPaCounterion : int
        total number of PA counterion
    numPcCounterion : int
        total number of PC counterion

    lenPa : int
        number of beads in a polyanion chain.
    lenPc : int
        number of beads in a polycation chain.
    numPaBeads : int
        total number of Beads in polyanion chains.
    numPa : int
        number of polyanion chain.
    numPc : int
        number of polycation chain.
    numChains : int
        number of polymer chains.
    numPcBeads : int
        total number of beads in polycation chains.
    numPaPcBeads : int
        total number of beads in PA and PC chains.

    chargeRatio : float
        charge ratio x.
    saltCon : float
        salt concentration.
    chargeFraction : float
        the fraction of beads bearing charges.
    numSalt : int
        number of add salt beads.

    numBeads : int
        Total number of beads including PA, PC, counterion and salt.
    
    !!!-----Considering using INHERITANCE --------!!!
    d : class dump object
        class object that holds dump file information.
    nsnaps : int
        number of snapshots in the whole d object.
    natoms : int
        number of atoms in one snapshots. should be consistent with setup up.

    xBox, yBox, zBox : float
        length of the Box in x,y,z direction.

    atomsInfo: 2d - ndarray
        Atoms information, Atoms ID, type and spatial information.
        sorted by Atoms's ID
    atomsId : 1d - ndarray, dtype = int64
        atoms ID information.
    atomsType : 1d - ndarray, dytpe = int64
        atoms Type information.
    atomsCoord : 2d - ndarray, dtype = float64
        atoms scaled coordination.

    molId : 1d - ndarray, dtype = int64
        molecule Id information
    numMols : int
        number of molecule number.

    clusterId : 1d - ndarray, dtype = int
        ID of cluster
    numClusters : int
        number of clusters
    clusterSize : 1d - ndarray, dtype = int
        number of beads in cluster not including couterions and salt.
    clusterMember : 2d - ndarray, dtype = int
        store information: chain members in a cluster, including itself.
    clusterChain : 1d - ndarray, dtype = int
        the number of chains in a cluster
    numNeighbors : 1d - ndarray, dtype = int
        number of neighbors.
    neighborListChain : 2d - ndarray, dtype = int
        store bonding information into neighbor_list
        this information is useful later on when moving chains
        not that the neighborList here is incomplete because we only
        perform cluster analysis when two chains are not found in the 
        same cluster already.
    counterionNeighbor : 1d - ndarray, dtype = int
        store the information of counterion's neighbor.

    
    """

    def __init__(self, ifile, lenPa, numPa, lenPc, numPc, saltCon):
        """
        Constructor.
        """
        self.chargeFraction = 0.5
        self.counterionSize = 0.5
        self.saltSize = 0.5
        self.beadSize = 1.0
        #self.distFactor = 2.0
        self.chargeRepeat = 2
        
        self.ifile = ifile
        self.lenPa = lenPa
        self.lenPc = lenPc
        self.numPa = numPa
        self.numPc = numPc
        #self.chargeRatio = chargeRatio
        self.saltCon = saltCon
        

        # derived member
        self.numPaBeads = self.lenPa * self.numPa
        self.numPcBeads = self.lenPc * self.numPc
        self.chargeRatio = float(self.numPcBeads)/float(self.numPaBeads)
        self.numPaPcBeads = self.numPaBeads + self.numPcBeads
        self.numChains = self.numPa + self.numPc 

        self.numPaCounterion = int(self.chargeFraction * self.numPaBeads)
        self.numPcCounterion = int(self.chargeFraction * self.numPcBeads)
        self.numCounterion = self.numPaCounterion + self.numPcCounterion
        self.numSalt = int(self.saltCon * self.numPaBeads)
        if self.numSalt%2 == 0:
            pass
        else:
            self.numSalt += 1

        self.numBeads = self.numPaPcBeads + self.numCounterion + self.numSalt

        self.lj_cutoff_case1 = self.counterionSize * np.power(2.0, 1.0/6.0)
        self.lj_cutoff_case2 = 0.5 * (self.beadSize + self.counterionSize) * np.power(2.0, 1.0/6.0)
        self.lj_cutoff_case3 = 0.5 * (self.beadSize + self.counterionSize) * np.power(2.0, 1.0/6.0)
        self.lj_cutoff_case4 = self.beadSize * np.power(2.0, 1.0/6.0)
        

    def readFile(self, index = 0):
        """
        Read dump file by dumpy.py.
        """
        self.d = dp.dump(self.ifile)
        self.nsnaps = self.d.nsnaps
        
        # choose a snapshot which will be analyzed here between 0 to nsnaps-1.
        # for now choose the first snapshots, e.g. index = 0!
        assert index <= self.nsnaps-1 and index >= 0
        self.snapsIndex = self.d.time()[index]
        self.natoms = self.d.snaps[index].natoms
        assert self.natoms == self.numBeads

        self.xBox = self.d.snaps[index].xhi - self.d.snaps[index].xlo
        self.yBox = self.d.snaps[index].yhi - self.d.snaps[index].ylo
        self.zBox = self.d.snaps[index].zhi - self.d.snaps[index].zlo
        self.xhi = self.d.snaps[index].xhi
        self.xlo = self.d.snaps[index].xlo
        self.yhi = self.d.snaps[index].yhi
        self.ylo = self.d.snaps[index].ylo
        self.zhi = self.d.snaps[index].zhi
        self.zlo = self.d.snaps[index].zlo
        self.box = np.array([self.xBox, self.yBox, self.zBox])

        
        #
        # lable Molecule ID number to each group, one PC/PA chain are one single molecule.
        # counterion and salt atom is considered as one molecule.
        #
        
        self.atomsInfo = sorted(self.d.snaps[index].atoms, key = lambda x:x[0])
        self.atomsType = np.array([self.atomsInfo[i][1] for i in range(0, self.natoms)])
        self.atomsType = self.atomsType.astype(int)
        self.atomsId = np.array([self.atomsInfo[i][0] for i in range(0, self.natoms)])
        self.atomsId = self.atomsId.astype(int)
        self.atomsCoord = np.array([self.atomsInfo[i][2:] for i in range(0, self.natoms)])
        self.molId = np.zeros(self.natoms, dtype = int)
        self.molId[:self.numPcBeads] = (self.atomsId[:self.numPcBeads] - 1)  / self.lenPc + 1
        self.molId[self.numPcBeads:self.numPaPcBeads] = self.numPc + (self.atomsId[self.numPcBeads:self.numPaPcBeads] - self.numPcBeads -1)/self.lenPa + 1
        self.molId[self.numPaPcBeads:] = self.atomsId[self.numPaPcBeads:] - self.numPaPcBeads + self.numPa + self.numPc

        self.numMols = self.natoms - self.numPaPcBeads + self.numPa + self.numPc

        assert self.numMols == self.molId[-1]

        self.atomsCharge = np.zeros(self.natoms, dtype = int)
        for i in range(self.natoms):
            if i < self.numPcBeads:
                if (i+1)%self.chargeRepeat != 0:
                    self.atomsCharge[i] = 0
                else:
                    self.atomsCharge[i] = 1
                
            elif i < self.numPaPcBeads:
                if (i+1)%self.chargeRepeat != 0:
                    self.atomsCharge[i] = 0
                else:
                    self.atomsCharge[i] = -1
            elif i < (self.numPaPcBeads + self.numPcCounterion):
                self.atomsCharge[i] =  -1
            elif i < (self.numPaPcBeads + self.numCounterion):
                self.atomsCharge[i] = 1
            elif i < (self.numPaPcBeads + self.numCounterion + self.numSalt/2):
                self.atomsCharge[i] = -1
            else:
                self.atomsCharge[i] = 1

        self.numBonds = self.numPc * (self.lenPc -1) + self.numPa * (self.lenPa - 1)
        self.bondList = []
        for i in range(self.numBonds):
            if i < self.numPc * (self.lenPc - 1):
                j = i/(self.lenPc - 1)
                k = i%(self.lenPc - 1)
                self.bondList.append([self.lenPc * j + k + 1, self.lenPc * j + k + 2])

            else:
                i = i - self.numPc * (self.lenPc - 1)
                j = i/(self.lenPa - 1)
                k = i%(self.lenPa - 1)                
                self.bondList.append([self.numPcBeads + self.lenPa*j + k + 1, self.numPcBeads + self.lenPa*j + k + 2])
            
                
                
                
        

    def analysis(self):
        """
        Cluster analysis for the system.
        initial cluster id of each molucule is the same as its molecular id.
        1st:
            We perform cluster analysis to PA and PC chains.
            We only perform it when two molecules are NOT in the same cluster already.
            merge two clusters if there are members of two clusters that are close enough to each other
            the cluster id is always the smallest molecular id of its members
            
        2nd:
            We perform cluster analysis to other counterions and salt.
        3rd:
            Reorder cluster Id.
        """
        
        
        self.clusterId = np.array(range(1, self.numMols+1))
        self.numNeighbors = np.zeros(self.numChains, dtype = int)
        self.neighborListChain = np.zeros((self.numChains, self.numChains), dtype = int)

        dist_factor = 2.0
        self.counterionNeighbor = np.zeros(self.numMols, dtype = int)
        for i in range(self.numChains-1):
            for j in range(i+1 , self.numChains):

                if self.clusterId[i] != self.clusterId[j]:
                   
                   
                    dij_min = self.distance(i, j)
                    dij_criterion = dist_factor * self.cutoff(i, j)

                    if dij_min <= dij_criterion:
                        
                        self.neighborListChain[i, self.numNeighbors[i]] = j
                        self.neighborListChain[j, self.numNeighbors[j]] = i
                        self.numNeighbors[i] += 1
                        self.numNeighbors[j] += 1

                        if self.clusterId[i] <= self.clusterId[j]:
                            cluster_temp = self.clusterId[j]
                            for m in range(0, self.numMols):
                                if self.clusterId[m] == cluster_temp:
                                    self.clusterId[m] = self.clusterId[i]
                        else:
                            cluster_temp = self.clusterId[i]
                            for m in range(0, self.numMols):
                                if self.clusterId[m] == cluster_temp:
                                    self.clusterId[m] = self.clusterId[j]

        # perform cluster analysis to seek for counterions that condense onto polymer chains
        # the factor that controls how close two mols needed to be in order to belong
        # to the same cluster
        dist_factor = 4.0
        for i in range(self.numChains):
            for j in range(self.numChains, self.numMols):
                if self.clusterId[i] != self.clusterId[j]:
                    
                    
                    dij_min = self.distance(i, j)
                    dij_criterion = dist_factor * self.cutoff(i, j)
                    if dij_min <= dij_criterion:
                        self.clusterId[j] = self.clusterId[i]
                        self.counterionNeighbor[j] = i
        
        # reorder cluster id from 1 to N, where N is the total number of clusters
        self.numClusters = 0
        self.clusterSize = np.zeros(self.numMols, dtype = int)
        self.clusterMember = np.zeros((self.numMols, self.numMols), dtype = int)
        self.clusterChain = np.zeros(self.numMols, dtype = int)
        for i in range(self.numChains):
            if self.clusterId[i] > self.numClusters:
                # find the starting chain member of a new cluster
                self.numClusters += 1
                #self.clusterSize[self.numClusters] = 0
                k = 0
                cluster_temp = self.clusterId[i]
                for m in range(self.numMols):
                    if self.clusterId[m] == cluster_temp:
                        self.clusterId[m] = self.numClusters
                        if m < self.numPc:
                            self.clusterSize[self.numClusters - 1] += self.lenPc
    
                            self.clusterMember[self.numClusters - 1,k] = m
                            k += 1
                        elif m < self.numChains:
                            self.clusterSize[self.numClusters - 1] += self.lenPa
                            
                            self.clusterMember[self.numClusters - 1,k] = m
                            k += 1
                        else:
                            self.clusterSize[self.numClusters - 1] += 1
                self.clusterChain[self.numClusters - 1] = k

        # if there are counterions, then put all the remaining counterions into the (N+1)-th cluster
        # note that some counterions might belong to the cluster formed by polymer chains
        # this phenomenon is called counterion condensation
        for i in range(self.numChains, self.numMols):
            if self.clusterId[i] > self.numClusters:
                self.clusterId[i] = self.numClusters + 1;
                self.clusterSize[self.numClusters] += 1

        assert sum(self.clusterSize) == self.numBeads
        
        
        # We call PA chain cluster as a strong cluster
        self.strongCluster = []
        for i in range(self.numClusters):
            #if self.clusterMember[i][0] <= self.numPa - 1:
            if any(list(i > self.numPc - 1 for i in self.clusterMember[i])):
                self.strongCluster.append(self.clusterMember[i][:self.clusterChain[i]])
        

    def moveCluster(self):
        """
        move all beads/counterions in the same cluster to the same spatial region
        without straddling over the boundary of the simulation box
        """
        self.centX = 0.5 * self.xBox
        self.centY = 0.5 * self.yBox
        self.centZ = 0.5 * self.zBox
        self.centBox = np.array([self.centX, self.centY, self.centZ])

        self.chainCenter = np.zeros(self.numChains,dtype = 'int')
        for m in range(self.numChains):
            dist_cent = float('inf')
            atom_start, atom_end = self.molindex2atomindex(m)
            for j in range(atom_start, atom_end +1):
                ds = self.atomsCoord[j] - self.centBox
                ds = ds - np.round(ds/self.box) * self.box
                rtocent = ds[0]*ds[0] + ds[1]*ds[1] + ds[2]*ds[2]
                if rtocent < dist_cent:
                    self.chainCenter[m] = j
                    dist_cent = rtocent


    
        self.atomsNewCoord = self.atomsCoord.copy()
        # move each chain to the same spatical region
        # use the bead closest to the box center to determin hwo the move
        # should be done.

        for m in range(self.numChains):
            j = self.chainCenter[m]
            self.map_back_bead(j)
            atom_start, atom_end = self.molindex2atomindex(m)
            assert j >= atom_start and j <= atom_end
            for i in range(j, atom_start, -1):
                self.move_bead(i-1,i)
            for i in range(j, atom_end):
                self.move_bead(i+1, i)

        # indentify the chain in a cluster that is closest to the center
        # of the simulation box

        self.clusterCenter = np.zeros(self.numClusters, dtype = 'int')
        for i in range(self.numClusters):
            dist_cent = float('inf')
            size_of_cluster = self.clusterChain[i]

            for k in range(size_of_cluster):
                m = self.clusterMember[i,k]
                assert self.clusterId[m] == i + 1 # Note that the cluster ID number is from 1 not 0!
                atom_start, atom_end = self.molindex2atomindex(m)
                chain_len = atom_end - atom_start + 1
                x_ave = sum(self.atomsNewCoord[atom_start:atom_end+1,0])/chain_len
                y_ave = sum(self.atomsNewCoord[atom_start:atom_end+1,1])/chain_len
                z_ave = sum(self.atomsNewCoord[atom_start:atom_end+1,2])/chain_len
                coord_ave = np.array([x_ave, y_ave, z_ave])
                ds = coord_ave - self.centBox
                ds = ds - np.round(ds/self.box) * self.box
                rtocent = ds[0]*ds[0] + ds[1]*ds[1] + ds[2]*ds[2]
                if rtocent < dist_cent:
                    self.clusterCenter[i] = m
                    dist_cent = rtocent


        # move chains belonging to the same cluster together in the same spatial region
        self.move_yesno = np.zeros(self.numChains, dtype = 'int')
        
        for i in range(self.numClusters):
            j = self.clusterCenter[i]
            self.map_back_chain(j)
            self.move_yesno[j] = 1
            size_of_cluster = self.clusterChain[i]

            if size_of_cluster > 1:
                move_cluster = 0
                while move_cluster == 0:
                    for k in range(size_of_cluster):
                        m = self.clusterMember[i,k]
                        assert self.clusterId[m] == i + 1

                        if self.move_yesno[m] == 0:
                            for n in range(self.numNeighbors[m]):
                                neighbor_chain_id = self.neighborListChain[m,n]
                                if self.move_yesno[neighbor_chain_id] == 1:
                                    self.move_mols(m, neighbor_chain_id)
                                    self.move_yesno[m] = 1
                                    break

                    move_cluster = 1
                    for k in range(size_of_cluster):
                        m = self.clusterMember[i,k]
                        assert self.clusterId[m] == i + 1
                        if self.move_yesno[m] == 0:
                            move_cluster = 0
                            break

        # move counterion to be close to its neighboring chains
        for i in range(self.numChains, self.numMols):
            if self.clusterId[i] <= self.numClusters:
                neighbor_chain_index = self.counterionNeighbor[i] # 
                self.move_mols(i, neighbor_chain_index)




                
# Start to write File: PDB file, dump file and a summary file



    def writeDumpFile(self):
        """
        Write a new dump file.
        """
        dirName = "./ClusterAnalysisPE.la{0}.na{1}.lc{2}.nc{3}.s{4}/".\
            format(self.lenPa, self.numPa, self.lenPc, self.numPc,\
                       self.saltCon)
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        dumpFileName = dirName + "New_dump_atom_pe_la{0}_lc{1}_npa{2}_cf{3}_s{4}_x{5}_snap{6}".\
                       format(self.lenPa, self.lenPc, self.numPa, self.chargeFraction,\
                              self.saltCon, self.chargeRatio, self.snapsIndex)
        ofile = open(dumpFileName, 'w')
        dumpFileString = "ITEM: TIMESTEP\n{0}\n".format(self.snapsIndex) + \
                         "ITEM: NUMBER OF ATOMS\n{0}\n".format(self.natoms) + \
                         "ITEM: BOX BOUNDS pp pp pp\n" + \
                         "{0} {1}\n".format(self.xlo, self.xhi) + \
                         "{0} {1}\n".format(self.ylo, self.yhi) + \
                         "{0} {1}\n".format(self.zlo, self.zhi) + \
                         "ITEM: ATOMS id type xs ys zs\n"
        ofile.write(dumpFileString)
        for i in range(self.natoms):
            scale = (self.atomsNewCoord[i] - np.array([self.xlo, self.ylo, self.zlo]) )/self.box
            content = "{0} {1} {2} {3} {4}\n".format(self.atomsId[i] , self.atomsType[i], \
                                                    scale[0], scale[1], scale[2])
            ofile.write(content)

        ofile.close()

    def writePDBFile(self):
        """
        write two PDB files.
        """
        dirName = "./ClusterAnalysisPE.la{0}.na{1}.lc{2}.nc{3}.s{4}/".\
            format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                       self.saltCon)
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        PDBFileName = dirName + "dump_atom_pe_la{0}_lc{1}_npa{2}_cf{3}_s{4}_x{5}_snap{6}.pdb".\
                       format(self.lenPa, self.lenPc, self.numPa, self.chargeFraction,\
                              self.saltCon, self.chargeRatio, self.snapsIndex)
        ofile = open(PDBFileName, 'w')
        PDBFileString = "CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f}{6}{7}\n".\
                        format(self.box[0], self.box[1], self.box[2], 90.0, 90.0, 90.0, \
                                " P 1       ", "   1" )
        ofile.write(PDBFileString)
        for m in range(self.numMols):
            atom_start, atom_end = self.molindex2atomindex(m)
            for i in range(atom_start, atom_end+1):
                aname = string.ascii_uppercase[self.atomsType[i] - 1]
                cid = self.clusterId[m]
                csize = np.round(self.clusterSize[cid-1]/100) # divide by 100 to conform with the format requirement of PDB file
                content  = "ATOM  {0:5}{1:>4}      {2:4}    {3:8.3f}{4:8.3f}{5:8.3f}  1.00{6:6.2f}\n".\
                          format(self.atomsId[i], aname, cid, self.atomsNewCoord[i,0], self.atomsNewCoord[i, 1],\
                                 self.atomsNewCoord[i, 2], csize)
                ofile.write(content)

        ofile.close()

        PDBFileName = dirName +  "dump_atom_pe_la{0}_lc{1}_npa{2}_cf{3}_s{4}_x{5}_snap{6}.nfc.pdb".\
                       format(self.lenPa, self.lenPc, self.numPa, self.chargeFraction,\
                              self.saltCon, self.chargeRatio, self.snapsIndex)
        ofile = open(PDBFileName, 'w')
        PDBFileString = "CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f}{6}{7}\n".\
                        format(self.box[0], self.box[1], self.box[2], 90.0, 90.0, 90.0, \
                                " P 1       ", "   1" )
        ofile.write(PDBFileString)
        
        for m in range(self.numMols):
            if self.clusterId[m] <= self.numClusters:
                atom_start, atom_end = self.molindex2atomindex(m)
                for i in range(atom_start, atom_end+1):
                    aname = string.ascii_uppercase[self.atomsType[i] - 1]
                    cid = self.clusterId[m]
                    csize = np.round(self.clusterSize[cid-1]/100) # divide by 100 to conform with the format requirement of PDB file
                    content  = "ATOM  {0:5}{1:>4}      {2:4}    {3:8.3f}{4:8.3f}{5:8.3f}  1.00{6:6.2f}\n".\
                               format(self.atomsId[i], aname, cid, self.atomsNewCoord[i,0], self.atomsNewCoord[i, 1],\
                                      self.atomsNewCoord[i, 2], csize)
                    ofile.write(content)
            
        ofile.close()    


    def writeRestartFile(self):
        """
        write restart data file
        """
        self.mass = []
        if int(self.saltCon) == 0:
            self.ntypes = 6
        else:
            self.ntypes = 8

        for i in range(self.ntypes):
            self.mass.append(1)

        self.bdtypes = 1
        self.angtypes = 0
        self.dihtypes = 0
        self.improtypes = 0
            
        dirName = "./ClusterAnalysisPE.la{0}.na{1}.lc{2}.nc{3}.s{4}/".\
                  format(self.lenPa, self.numPa, self.lenPc, self.numPc,\
                             self.saltCon)
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        iFileLammpsName = dirName +  "res.data.pe.la{0}.na{1}.lc{2}.nc{3}.r{4}.s{5}.lammps".\
                          format(self.lenPa, self.numPa, self.lenPc, self.numPc, self.chargeRepeat, self.saltCon)
        iFileLammps = open(iFileLammpsName, 'w')

        iFileLammpsHeader = "restart data file for mixtures for charged polymer chains\n" +\
                            "\n" + \
                            "{0:10d}    atoms\n".format(self.natoms) + \
                            "{0:10d}    bonds\n".format(self.numBonds) + \
                            "{0:10d}    angles\n".format(0) + \
                            "{0:10d}    dihedrals\n".format(0) + \
                            "{0:10d}    impropers\n".format(0) + \
                            "\n" +\
                            "{0:10d}    atom types\n".format(self.ntypes) + \
                            "{0:10d}    bond types\n".format(self.bdtypes) + \
                            "{0:10d}    angle types\n".format(self.angtypes) + \
                            "{0:10d}    dihedral types\n".format(self.dihtypes) + \
                            "{0:10d}    improper types\n".format(self.improtypes) + \
                            "\n" + \
                            " {0:16.8f} {1:16.8f}   xlo xhi\n".format(self.xlo, self.xhi) + \
                            " {0:16.8f} {1:16.8f}   ylo yhi\n".format(self.ylo, self.yhi) + \
                            " {0:16.8f} {1:16.8f}   zlo zhi\n".format(self.zlo, self.zhi) + \
                            "\n" + \
                            "Masses\n" + \
                            "\n"

        iFileLammps.write(iFileLammpsHeader)
        for i in range(self.ntypes):
            iFileLammps.write( "{0}  {1:8.3f}\n".format(i+1, self.mass[i]))

        iFileLammps.write("\nAtoms\n\n")

        for i in range(self.natoms):
            iFileLammps.write("{0} {1} {2} {3} {4} {5} {6}\n".format(i+1, \
                                                                   self.molId[i], \
                                                                   self.atomsType[i], \
                                                                   self.atomsCharge[i], \
                                                                   self.atomsNewCoord[i][0], \
                                                                   self.atomsNewCoord[i][1], \
                                                                   self.atomsNewCoord[i][2]))

        iFileLammps.write("\nBonds\n\n")
        for i in range(self.numBonds):
            iFileLammps.write("{0} 1 {1} {2}\n".format(i+1, self.bondList[i][0], self.bondList[i][1]))

        
        iFileLammps.close()

        
        

    def writeSummaryFile(self):
        dirName = "./ClusterAnalysisPE.la{0}.na{1}.lc{2}.nc{3}.s{4}/".\
            format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                       self.saltCon)       
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        summaryFileName =  dirName + "dump_atom_pe_la{0}_lc{1}_npa{2}_cf{3}_s{4}_x{5}_snap{6}.cluster".\
                       format(self.lenPa, self.lenPc, self.numPa, self.chargeFraction,\
                              self.saltCon, self.chargeRatio, self.snapsIndex)
        ofile = open(summaryFileName, 'w')
        ####################################################################################
        # ofile.write("# a single PA or PC chain can be defined as cluster!\n")            #
        # ofile.write("Number of Clusters: " + str(self.numClusters) + '\n')               #
        #                                                                                  #
        # for i in range(self.numClusters):                                                #
        #     ofile.write(str(self.clusterChain[i]) + '\n')                                #
        #                                                                                  #
        # ofile.write('\n')                                                                #
        # ofile.write('# Only a PA chain can be defined as a strong cluster!\n' + \        #
        #             '# Thus it can have at most numPa clusters in the simulation box\n') #
        # ofile.write("Number of strong clusters: " + str(len(self.strongCluster)) + '\n') #
        # for i in self.strongCluster:                                                     #
        #     for j in i:                                                                  #
        #         ofile.write(str(j) + ' ')                                                #
        #     ofile.write('\n')                                                            #
        #                                                                                  #
        # ofile.write('\n' + '\n Done!')                                                   #
        ####################################################################################
        headline = "Summary of Cluster Analysis:\n" +\
                   "\n" +\
                   "System Configuration:\n" +\
                   "Polycation(PC) chain length: {0}\n".format(self.lenPc) +\
                   "Numbers of polycations: {0}\n".format(self.numPc) +\
                   "Polyanion(PA) chain length: {0}\n".format(self.lenPa) +\
                   "Numbers of polyanion: {0}\n".format(self.numPa) +\
                   "Charge ratio: {0}\n".format(self.chargeRatio) +\
                   "Salt concentration/Polyanion monomer concentration: {0}\n".format(self.saltCon) +\
                   "Number of total beads: {0}\n".format(self.numBeads) +\
                   "monomer density: {0}\n\n".format(self.numPaPcBeads/self.xBox/self.yBox/self.zBox)
        ofile.write(headline)
        notice = "- - - - - - - - - - - - - - - - \n" +\
                 "Chain ID info: chain's ID which is less than {0} are polycation chain.\n".format(self.numPc) +\
                 "               chain's ID which is equal to or greater than {0} are polyanion chain\n.".format(self.numPc) +\
                 "Definition about different phases:\n" +\
                 "\phi_1 is a single polyanion chain or a single polycation chain, AKA soluble phase\n" +\
                 "\phi_2 is a small cluster with only one PA chain and several PC chains\n" +\
                 "\phi_3 is a large cluster with multiple condensed PA and PC chains\n" +\
                 "- - - - - - - - - - - - - - - - \n"
        ofile.write(notice)
        
        # phi1 is the a single PA or PC chain, as known as soluble phase
        # phi2 is a small cluster with only one PA chain and several PC chains
        # phi3 is a large cluster with multiple PA anc PC chains
        phi1, phi2, phi3 = 0, 0, 0

        for i in range(self.numClusters):
            if self.clusterChain[i] == 1:
                phi1 += 1
                ofile.write("ph1: {0}\n").format(model.clusterMember[i][0])
            else:
                pcchains = 0
                for j in range(self.clusterChain[i]):
                    if self.clusterMember[i][j] >= self.numPc:
                        pcchains += 1
                if pcchains >= 2:
                    phi3 += 1
                    ofile.write("phi3: ")
                    for j in range(self.clusterChain[i]):
                        ofile.write(str(self.clusterMember[i][j]) + " ")
                    ofile.write("\n")
                else:
                    phi2 += 1
                    ofile.write("phi2: ")
                    for j in range(self.clusterChain[i]):
                        ofile.write(str(self.clusterMember[i][j]) + " ")
                    ofile.write("\n")

        last = "- - - - - - - - - - - - - - - - \n" +\
               "Number of different phases in the simulation: \n" +\
               "\phi_1: {0}\n".format(phi1) +\
               "\phi_2: {0}\n".format(phi2) +\
               "\phi_3: {0}\n".format(phi3) +\
               "Total phases: {0}\n".format(phi1 + phi2 + phi3) +\
               "Total clusters: {0}\n".format(self.numClusters)
        ofile.write(last)
        ofile.close()           # 

        
        
# -------This part needs to be fixed ----
# Since salt atom size is not defined and distance factor is also unclear.
                    
    def cutoff(self, a, b):
        """
        the cutoff bwtween two molecules used to judge if they are clustered togother
        or not       
        """
        assert a >= 0 and a < self.numMols
        assert b >= 0 and b < self.numMols

        if a > self.numChains - 1:
            if b > self.numChains - 1:
                lj_cutoff = self.lj_cutoff_case1
            else:
                lj_cutoff = self.lj_cutoff_case2
        else:
            if b > self.numChains - 1:
                lj_cutoff = self.lj_cutoff_case3
            else:
                lj_cutoff = self.lj_cutoff_case4

        return lj_cutoff

    
    def distance(self, a, b):
        """
        distance between two molecules defined asd the minimum 
        separation of any two atoms respectively from the two molecules.
        """
        a_start, a_end = self.molindex2atomindex(a)
        b_start, b_end = self.molindex2atomindex(b)

        # A smart way considering periodic boundary condition
        # to compute distance of two atoms by Shengfeng Cheng
        
        r_min = float('inf')
        for ia in range(a_start, a_end + 1):
            for ib in range(b_start, b_end + 1):
                ds = self.atomsCoord[ia] - self.atomsCoord[ib]
                ds = ds - np.round(ds/self.box) * self.box
                r_temp = ds[0] * ds[0] + ds[1] * ds[1] + ds[2] * ds[2]
                if r_temp < r_min:
                    r_min = r_temp

        return np.sqrt(r_min)
    

    def molindex2atomindex(self, x):
        """
        Translate the molecule index to Atom index.
        Return a tuple of begin atom index and end of atom index.
        """
        assert x >= 0 and x < self.numMols

        if x < self.numPc:
            start = x * self.lenPc
            end = (x + 1) * self.lenPc - 1 
        elif x < self.numChains:
            start = self.numPcBeads + (x - self.numPc) * self.lenPa
            end = self.numPcBeads + (x - self.numPc + 1) * self.lenPa - 1
        else:
            start = self.numPaPcBeads + x - self.numChains
            end = start

        return start, end


    def map_back_bead(self, a):
        """
        map a bead back to the simulation box.
        """
        ds = self.atomsNewCoord[a] - self.centBox
        ds = ds - np.round(ds/self.box) * self.box
        self.atomsNewCoord[a] = self.centBox + ds

    def move_bead(self, a, b):
        """
        move two beads to the same spatial region.
        """
        ds = self.atomsNewCoord[a] - self.atomsNewCoord[b]
        ds = ds - np.round(ds/self.box) * self.box
        self.atomsNewCoord[a] = self.atomsNewCoord[b] + ds


    def map_back_chain(self, a):
        """
        map a chain back to  the simulation box.
        """
        a_start, a_end = self.molindex2atomindex(a)
        chain_len_temp = a_end - a_start + 1
        x_ave_temp = sum(self.atomsNewCoord[a_start:a_end+1,0])/chain_len_temp
        y_ave_temp = sum(self.atomsNewCoord[a_start:a_end+1,1])/chain_len_temp
        z_ave_temp = sum(self.atomsNewCoord[a_start:a_end+1,2])/chain_len_temp

        ave_temp = np.array([x_ave_temp, y_ave_temp, z_ave_temp])
        ds = ave_temp - self.centBox
        ds = ds - np.round(ds/self.box) * self.box
        self.atomsNewCoord[a_start:a_end+1] += self.centBox - ave_temp + ds

    def move_mols(self, a, b):
        """
        map a molecule to the target molecule.
        """
        a_start, a_end = self.molindex2atomindex(a)
        b_start, b_end = self.molindex2atomindex(b)

        chain_len_temp = b_end - b_start + 1
        x_target = sum(self.atomsNewCoord[b_start:b_end +1,0])/chain_len_temp
        y_target = sum(self.atomsNewCoord[b_start:b_end +1,1])/chain_len_temp
        z_target = sum(self.atomsNewCoord[b_start:b_end +1,2])/chain_len_temp

        chain_len_temp = a_end - a_start + 1
        x_ave_temp = sum(self.atomsNewCoord[a_start:a_end +1,0])/chain_len_temp
        y_ave_temp = sum(self.atomsNewCoord[a_start:a_end +1,1])/chain_len_temp
        z_ave_temp = sum(self.atomsNewCoord[a_start:a_end +1,2])/chain_len_temp

        target = np.array([x_target, y_target, z_target])
        ave_temp = np.array([x_ave_temp, y_ave_temp, z_ave_temp])

        ds = ave_temp - target
        ds = ds - np.round(ds/self.box) * self.box

        self.atomsNewCoord[a_start:a_end+1] += - ave_temp + target + ds

    
        
        

if __name__ == "__main__":
    dirName = "./ClusterAnalysisPE.la{0}.na{1}.lc{2}.nc{3}.s{4}/".\
        format(len_pan, num_pan, len_pcat, num_pcat,\
                   salt_concentration)
    if not os.path.exists(dirName):
        os.makedirs(dirName)

    Sum = []
    SumFileName = dirName + sys.argv[1][-10:] +  'summary.txt'
    ofile = open(SumFileName, 'w')
    
    for i in range(1):
        
    
        model = ClusterAnalysis(filename, len_pan, num_pan, len_pcat, num_pcat, salt_concentration)
        model.readFile(i)
        start = time.time()
        model.analysis()
        model.moveCluster()
        model.writeDumpFile()
        model.writePDBFile()
        model.writeRestartFile()
        model.writeSummaryFile()
        Sum.append(model.numClusters)
        ofile.write(str(model.numClusters) + '\n')
        print time.time() - start

    Sum = np.array(Sum)
    SumMeanAndStd = str(Sum.mean()) + '\t' +  str(Sum.std()) + '\n'
    ofile.write(SumMeanAndStd)
    ofile.close()


