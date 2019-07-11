#!/usr/bin/python
# PolyelectrolyteBuilder.py
# Author: Yanfei Tang & Dr. Shengfeng Cheng

oneline = "This is a script that generates polyelectrolyte solution configuration for molecular dynamics simulation."

docstr = """
This is the script that generates the starting configuration of polyelectrolytes solution for molecular dynamics
simulation.
Essentially, binary polyelectrolytes solution contains polyanion chain and polycation chain may or may not
contain salt. The dynamics of such solution is not yet understood. The MD simulation uses bead-spring model to 
study such solution. For example, the sodium polystyrene sulfonate (PSSNa) is a polyanion, sodium ion is considered 
a positively charged couterion, sulfonate is considered a negatively charged beads. 

Usage:
./PolyelectrolyteBuilder.py len_pan num_pan len_pcat num_pcat vol_ratio charge_fraction salt_concentration
"""

# -----Library import -----
import numpy as np
import sys
import time
import string
import os

# -----------------------


# class definition

class PEBuilder:
    """
    The Attribute is explained in every function when it firstly apperears.
    """

    def __init__(self, len_pan, num_pan, len_pcat, num_pcat, vol_ratio, charge_fraction, salt_concentration):
        """
        Constructor.
        Attributes:
        -----------
        lenPa: int
            length of Polyanion, eg 128 monomer/ 128 beads per polyanion chain in bead-spring model
        numPa: int
            number of polyanion chains.
        lenPc: int
            lenght of polycation, eg 32 monomer/ 32 beads in one polycation chain.
        numPc: int
            number of polycation chains.
        volRatio: float
            staring volume ration of polymer chains, was 0.03 in Dobrynin's 2004 PRL paper. 
            Molecular Dynamics Simulations of Electrostatics Layer-by_Layer Self-Assembly
        chargeFraction: float
            the fraction of beads bearing charges
        saltCon: float
            saltCon is actually a ratio (not confused by it name concentration)
            A ration btw the number of add positively ions (monovalently charged for now) and
            the number of negatively charged (-e) polymeric subunits on polyanion chains.
        chargeRepeat: int
            the reciprocal of the chargeFraction
        totVol: float
            the total volume of the simulation box.
        xBox, yBox, zBox, box: float, np.array
            the length of the simulation box.
        segment: float
            segment length on each polymer chain, in the unit of \sigma
            this length is for standard FENE potential
        stiff: float
            parameter controls the stiffness of chains, distance btw No.1 and No.3 in 1-2-3 configuration
            = 2 * segment for rigid straight chains
        lx, ly, lz, hx, hy, hz: float
            self-explained
        """
        self.lenPa = len_pan
        self.numPa = num_pan
        self.lenPc = len_pcat
        self.numPc = num_pcat
        self.volRatio = vol_ratio
        self.chargeFraction = charge_fraction
        self.saltCon = salt_concentration
        self.chargeRepeat = int(np.floor(1/self.chargeFraction))
        self.totVol = (self.numPa * self.lenPa + self.numPc * self.lenPc)/self.volRatio
        self.xBox = np.power(self.totVol, 1/3.0)
        self.yBox = self.xBox
        self.zBox = self.xBox
        self.box = np.array([self.xBox, self.yBox, self.zBox])
        self.segment = 0.9609
        self.stiff = 1.4 * self.segment
        self.stiffangle = 2 * np.arcsin(0.5 * self.stiff/self.segment)

        self.lx, self.ly, self.lz = 0.0, 0.0, 0.0
        self.lxyz = np.array([self.lx, self.ly, self.lz])
        self.hx, self.hy, self.hz = self.lx + self.xBox, self.ly + self.yBox, self.lz + self.zBox
        self.hxyz = np.array([self.hx, self.hy, self.hz])
        self.dstot = np.zeros(3, dtype = 'float')


    def genChains(self):
        """
        This is a routine to generate chains.
        Attributes:
        -----------
        numMonomer: int
            number of monomers, thus number of beads in chains
        numBonds: int
            number of bonds
        numMols: int
            number of molecules
        numCations: int
            nuber of cations
        numAnions: int
            number of anions
        beadCharge: int
            +1 or -1 charge
        beadType: int
            type of the bead
        atomsCoords: list
            the coordinations of atoms
        atomsType: list
            the list contains the atoms' type
        atomsCharge: list
            the list contains the atoms' charge information
        molId: list
            the list contains the atoms' ID information
        """
        self.numMonomer = 0
        self.numBonds = 0
        self.numMols = 0
        self.numCations = 0
        self.numAnions = 0

        self.atomsCoords = []
        self.atomsType = []
        self.atomsCharge = []
        self.molId = []
        self.bondList = []
        
        for i in range(self.numPa + self.numPc):

            if i < self.numPc:
                # polycation chains, charge in LJ units of LAMMPS
                # electron charge would be 10.54 using bare LAMMPS LJ units
                # the dielectric constans of solvent is effectively taken as 111 when assign 1 to +e
                # just need to set dielectric as 0.72 in LAMMPS ot mimic water with dielectric constant 80
                self.beadCharge = 1
                self.beadType = 1 # atomic type for neutral beads in polycation chains
                self.chain = self.lenPc
            else:
                self.beadCharge = -1 # polyanion chains
                self.beadType = 3 # atomic type for neutral beads in polyanion chains
                self.chain = self.lenPa

            self.numMols += 1

            # generate the first bead of each chain randomly
            self.numMonomer += 1
            self.cxyz = np.random.rand(3) * self.box + self.lxyz

            self.atomsCoords.append(self.cxyz)
            #self.atomsType.append(self.beadType)

            # decide if the first bead is charged or not
            if self.chargeRepeat == 1:
                self.atomsCharge.append(self.beadCharge)
                self.atomsType.append(self.beadType + 1)
                if i < self.numPc:
                    self.numCations += 1
                else:
                    self.numAnions += 1
            else:
                self.atomsType.append(self.beadType)
                self.atomsCharge.append(0)

            self.molId.append(self.numMols)

            self.currpxyz = self.cxyz

            # follow random walk to generate the chain
            # generate the seconb bead of the chain
            self.theta, self.phi = np.random.rand(2) * np.array([np.pi, 2 * np.pi])
            self.ds = np.array([np.cos(self.theta), np.sin(self.theta) * np.cos(self.phi),\
                                np.sin(self.theta) * np.sin(self.phi)]) * self.segment

            self.cxyz = self.currpxyz + self.ds

            self.numMonomer += 1
            self.atomsCoords.append(self.cxyz)

            # decide if the second bead is charged or not
            if 2%self.chargeRepeat == 0:
                self.atomsCharge.append(self.beadCharge)
                self.atomsType.append(self.beadType + 1)
                if i < self.numPc:
                    self.numCations += 1
                else:
                    self.numAnions += 1
            else:
                self.atomsCharge.append(0)
                self.atomsType.append(self.beadType)

            self.molId.append(self.numMols)
            
            self.numBonds += 1
            self.bondList.append([self.numMonomer - 1, self.numMonomer])

            self.currpxyz = self.cxyz

            self.currtheta = self.theta
            self.currphi = self.phi

            self.dstot += self.ds

            # generating the rest beads of the chain

            for k in range(3, self.chain+1):
                # only accept atoms that are beyong certain distance
                # from the atom precding the current atom in the chain
                self.theta, self.phi = np.random.rand() * np.array([np.pi - self.stiffangle, \
                                                                    2 * np.pi])
                self.ds1 = np.array([np.cos(self.theta), np.sin(self.theta) * np.cos(self.phi),\
                                np.sin(self.theta) * np.sin(self.phi)]) * self.segment

                self.reverseXZrotation()
                self.cxyz = self.currpxyz + self.ds

                self.numMonomer += 1
                self.atomsCoords.append(self.cxyz)

                if k % self.chargeRepeat == 0:
                    self.atomsCharge.append(self.beadCharge)
                    self.atomsType.append(self.beadType + 1)
                    if i < self.numPc:
                        self.numCations += 1
                    else:
                        self.numAnions += 1
                else:
                    self.atomsCharge.append(0)
                    self.atomsType.append(self.beadType)

                self.molId.append(self.numMols)
                self.numBonds += 1
                self.bondList.append([self.numMonomer - 1, self.numMonomer])

                self.currpxyz = self.cxyz

                self.currtheta = np.arccos(self.ds[0]/self.segment)
                if self.ds[2] > 0:
                    self.currphi = np.arccos(self.ds[1]/self.segment/np.sin(self.currtheta))
                else:
                    self.currphi = 2*np.pi - np.arccos(self.ds[1]/self.segment/np.sin(self.currtheta))

                self.dstot += self.ds

        print "%d beads are generated.\n" % self.numMonomer 
        assert self.numMonomer == self.numPc * self.lenPc + self.numPa * self.lenPa, \
            "The number of monomers in chains is wrong!\n"
        assert self.numCations == int(np.floor(self.lenPc * self.chargeFraction)*self.numPc), \
        "The number of positively charged beads is wrong!\n"
        assert self.numAnions == int(np.floor(self.lenPa * self.chargeFraction)*self.numPa), \
        "The number of negatively charged beads is wrong!\n"


    def genCounterions(self):
        """
        start generating counterions.
        Attribute:
        ----------
        numCounterions: int
            number of counterions.
        """
        self.numCounterions = self.numCations + self.numAnions
        for i in range(self.numCounterions):
            self.numMols += 1
            self.cxyz = self.lxyz + self.box * np.random.rand(3)

            self.atomsCoords.append(self.cxyz)
            self.molId.append(self.numMols)

            if i < self.numCations:
                self.atomsCharge.append(-1) # negatively charged counterions for polycation chains
                self.atomsType.append(5)
            else:
                self.atomsCharge.append(1) # positively charged counterions for polyanion chains
                self.atomsType.append(6)

        print "%d counterions (including positive and negative) are generated.\n" % self.numCounterions


    def genSalt(self):
        """
        start generating salt ions.
        Attributes:
        -----------
        numSalt: int
            number of salt atoms
        """
        self.numSalt = int(2 * self.numAnions * self.saltCon) # the factor 2 counts for positive and negative ions from the added salt, and main electronutrality.
        if self.numSalt%2 == 0:
            pass
        else:
            self.numSalt += 1

        for i in range(self.numSalt):
            self.numMols += 1
            self.cxyz = self.lxyz + self.box * np.random.rand(3)

            self.atomsCoords.append(self.cxyz)
            self.molId.append(self.numMols)

            if i < self.numSalt/2:
                self.atomsCharge.append(-1) # negative ions from the added salt (monovalent)
                self.atomsType.append(7)
            else:
                self.atomsCharge.append(1) # positive ions from the added salt
                self.atomsType.append(8)
            
        print "%d salt ions (including positve and negative) are generated.\n" % self.numSalt


    def writeFiles(self, directory = "./"):
        """
        write data file and xyz file.
        Attribute:
        ----------
        """
        self.mass = []
        self.zero = 0
        self.natoms = self.numMonomer + self.numCounterions + self.numSalt
        self.nangles = 0
        self.ndihedrals = 0

        if int(np.ceil(self.saltCon)) == 0:
            self.ntypes = 6
        else:
            self.ntypes = 8

        # set masses of all beads to be 1
        # in principle, the mass of counterions and salt ions should be smaller
        # expect this difference will no matter in terms of complexation of polyelectrolytes
        for i in range(self.ntypes):
            self.mass.append(1)



        self.bdtypes = 1
        self.angtypes = 0
        self.dihtypes = 0
        self.improtypes = 0

        iFileLammpsName = directory + "data.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.lammps".\
                          format(self.lenPa, self.numPa, self.lenPc, self.numPc, self.volRatio, self.chargeRepeat, self.saltCon)
        iFileLammps = open(iFileLammpsName, 'w')

        iFileXYZName = directory + "data.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.xyz".\
                       format(self.lenPa, self.numPa, self.lenPc, self.numPc, self.volRatio, self.chargeRepeat, self.saltCon)
        iFileXYZ = open(iFileXYZName, 'w' )

        iFileXYZ.write("{0}\n".format(self.natoms))
        iFileXYZ.write("data.polyelectrolyte.xyz\n")

        iFileLammpsHeader = "data file for mixtures of charged polymer chains\n" + \
                            "\n" + \
                            "{0:10d}    atoms\n".format(self.natoms) + \
                            "{0:10d}    bonds\n".format(self.numBonds) + \
                            "{0:10d}    angles\n".format(self.nangles) + \
                            "{0:10d}    dihedrals\n".format(self.ndihedrals) + \
                            "{0:10d}    impropers\n".format(self.zero) + \
                            "\n" +\
                            "{0:10d}    atom types\n".format(self.ntypes) + \
                            "{0:10d}    bond types\n".format(self.bdtypes) + \
                            "{0:10d}    angle types\n".format(self.angtypes) + \
                            "{0:10d}    dihedral types\n".format(self.dihtypes) + \
                            "{0:10d}    improper types\n".format(self.improtypes) + \
                            "\n" + \
                            " {0:16.8f} {1:16.8f}   xlo xhi\n".format(self.lx, self.hx) + \
                            " {0:16.8f} {1:16.8f}   ylo yhi\n".format(self.ly, self.hy) + \
                            " {0:16.8f} {1:16.8f}   zlo zhi\n".format(self.lz, self.hz) + \
                            "\n" + \
                            "Masses\n" + \
                            "\n"

        iFileLammps.write(iFileLammpsHeader)
        for i in range(self.ntypes):
            iFileLammps.write( "{0}  {1:8.3f}\n".format(i+1, self.mass[i]))

        iFileLammps.write("\nAtoms\n\n")
        
        

        for i in range(self.natoms):
            if self.atomsType[i] == 1 or self.atomsType[i] == 3:
                iFileXYZ.write("S {0} {1} {2}\n".format(self.atomsCoords[i][0], \
                                                        self.atomsCoords[i][1], \
                                                        self.atomsCoords[i][2]))
            elif self.atomsType[i] == 2:
                iFileXYZ.write("P {0} {1} {2}\n".format(self.atomsCoords[i][0], \
                                                        self.atomsCoords[i][1], \
                                                        self.atomsCoords[i][2]))
            elif self.atomsType[i] == 4:
                iFileXYZ.write("N {0} {1} {2}\n".format(self.atomsCoords[i][0], \
                                                        self.atomsCoords[i][1], \
                                                        self.atomsCoords[i][2]))
            elif self.atomsType[i] == 5:
                iFileXYZ.write("A {0} {1} {2}\n".format(self.atomsCoords[i][0], \
                                                        self.atomsCoords[i][1], \
                                                        self.atomsCoords[i][2]))
            elif self.atomsType[i] == 6:
                iFileXYZ.write("C {0} {1} {2}\n".format(self.atomsCoords[i][0], \
                                                        self.atomsCoords[i][1], \
                                                        self.atomsCoords[i][2]))
            elif self.atomsType[i] == 7:
                iFileXYZ.write("I {0} {1} {2}\n".format(self.atomsCoords[i][0], \
                                                        self.atomsCoords[i][1], \
                                                        self.atomsCoords[i][2]))
            elif self.atomsType[i] == 8:
                iFileXYZ.write("K {0} {1} {2}\n".format(self.atomsCoords[i][0], \
                                                        self.atomsCoords[i][1], \
                                                        self.atomsCoords[i][2]))

            iFileLammps.write("{0} {1} {2} {3} {4} {5} {6}\n".format(i+1, \
                                                                   self.molId[i], \
                                                                   self.atomsType[i], \
                                                                   self.atomsCharge[i], \
                                                                   self.atomsCoords[i][0], \
                                                                   self.atomsCoords[i][1], \
                                                                   self.atomsCoords[i][2]))

        iFileLammps.write("\nBonds\n\n")
        for i in range(self.numBonds):
            iFileLammps.write("{0} 1 {1} {2}\n".format(i+1, self.bondList[i][0], self.bondList[i][1]))

        iFileXYZ.close()
        iFileLammps.close()
                
        

        
            

# ---- auxillary function -------
    def reverseXZrotation(self):
        """
        Maps the displacement vector in the current frame
        back to the original frame of the simulation box.
        """
        rot = np.zeros((3, 3), dtype = 'float64')
        rot[0, 0] = np.cos(self.currtheta)
        rot[0, 1] = -np.sin(self.currtheta)
        rot[0, 2] = 0.0
        rot[1, 0] = np.sin(self.currtheta) * np.cos(self.currphi)
        rot[1, 1] = np.cos(self.currtheta) * np.cos(self.currphi)
        rot[1, 2] = - np.sin(self.currphi)
        rot[2, 0] = np.sin(self.currtheta) * np.sin(self.currphi)
        rot[2, 1] = np.cos(self.currtheta) * np.sin(self.currphi)
        rot[2, 2] = np.cos(self.currphi)

        self.ds  = np.dot(rot, self.ds1)

        


if __name__ == "__main__":
    len_pan = int(sys.argv[1])
    num_pan = int(sys.argv[2])
    len_pcat = int(sys.argv[3])
    num_pcat = int(sys.argv[4])
    vol_ratio = float(sys.argv[5])
    charge_fraction = float(sys.argv[6])
    salt_concentration = float(sys.argv[7])
    model = PEBuilder(len_pan, num_pan, len_pcat, num_pcat, vol_ratio, charge_fraction,\
                      salt_concentration)
    model.genChains()
    model.genCounterions()
    model.genSalt()
    model.writeFiles()



    
