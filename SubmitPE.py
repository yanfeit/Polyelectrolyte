#!/usr/bin/python
# submitPE.py

oneline = "This is a script that used to submit job on tbird."

docstr = """
Binary polyelectrolytes system in solution.
Author: Yanfei Tang @ VT Physics Dept.
Date: May 26, 2016


Job summission procedure:
- PBS bash file:
--- walltime, nodes, total cores, excutable lammps name
--- input file, log file, outfile


- lammps input file:
-- in_pe_pushoff
-- in_pe_equi

- data.file res.file.

- PolyelectrolytsBuilder.py

To Do List:


"""

import subprocess
import os
import platform
import time
import sys
import string
import numpy as np
import PolyelectrolyteBuilder as pe


class Job():
    """
    class Job means a single run on tbird. 
    """

    def __init__(self, walltime = 240,
                 nodes = 1,
                 ppn = 16,
                 lammps = "~/bin/lmp_tbird_15MAY2015",
                 len_pan = 128,
                 num_pan = 32,
                 len_pcat = 32,
                 num_pcat = 128,
                 vol_ratio = 0.01,
                 charge_fraction = 0.5,
                 salt_concentration = 1.0,
                 pretimestep = 0,
                 timestep = 2000000,
                 pushoff = True,
                 first = True
                 ):
        """
        A single job is only for a specific binary polyelectrolytes system.
        We need some parameter to specify the system, e.g. number of polyanion
        in the system, number of polycation in the system and so on. Also, we
        will need all the running are all put in one directory.

        The constructor is basicly self-explained.But...
        walltime: the maximum time allowed on PBS file
        nodes: number of nodes used.
        lammps: the excutable lammps name
        lenPa: length of polyanion
        numPa: number of polyanion
        lenPc: lenght of polycation
        numPc: number of polycation
        volRatio: essentially it is monomer density
        chargeFraction: how often the beads are charged
        saltCon: salt concentration
        timestep: the expected timestep for this job. It means after
        the job is finished, this valued is the last timestep.
        pushoff: is this job for pushoff purpose?
        first: is this job for the first time equilibrium?
        It is useful to write the input file.
        """
        self.walltime = walltime
        self.nodes = nodes
        self.ppn = ppn
        self.cores = nodes*ppn
        self.lammps = lammps
        self.lenPa = len_pan
        self.numPa = num_pan
        self.lenPc = len_pcat
        self.numPc = num_pcat
        self.volRatio = vol_ratio
        self.chargeFraction = charge_fraction
        self.saltCon = salt_concentration
        self.chargeRepeat = int(np.floor(1/self.chargeFraction))
        self.timestep = timestep
        self.pretimestep = pretimestep
        self.pushoff = pushoff
        self.first = first
        self.directory = "la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}/".\
                         format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                self.volRatio, self.chargeRepeat, self.saltCon)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        
        if self.pushoff:
            self.writeDataFile()

        self.writePushoffInputFile()

        if not self.pushoff:
            self.writeEquiInputFile()

        self.writePBSFile()

        print "Job info:\n"
        print "Length of polyanion: {0}\n".format(self.lenPa)
        print "Number of polyanion: {0}\n".format(self.numPa)
        print "Lenght of polycation: {0}\n".format(self.lenPc)
        print "Number of polycation: {0}\n".format(self.numPc)
        print "Monomer Density: {0}\n".format(self.volRatio)
        print "Every {0} bead is charged on polyelectrolyte\n".format(self.chargeRepeat)
        print "Salt concentration: {0}\n".format(self.saltCon)
        

        time.sleep(5)

        


    def writeDataFile(self):
        """
        We use PolyelectrolyteBuilder.py to create a initial system.
        Surely, we will need a pushoff to redistribute the system to 
        make it physical.
        """
        model = pe.PEBuilder(self.lenPa, self.numPa, self.lenPc, self.numPc,\
                             self.volRatio, self.chargeFraction, self.saltCon)
        model.genChains()
        model.genCounterions()
        model.genSalt()
        model.writeFiles(self.directory)
        self.dataFileName = "data.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.lammps".\
                            format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                   self.volRatio, self.chargeRepeat, self.saltCon)

    def writePushoffInputFile(self):
        """
        write a in_pe_pushoff input file.
        """
        self.in_pe_pushoffName = "la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.pushoff.txt".\
                                 format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                        self.volRatio, self.chargeRepeat, self.saltCon)
        
        self.in_pe_pushoffPath = self.directory + self.in_pe_pushoffName
        
        self.pushoffDumpFileName = "dump.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.pushoff".\
                                   format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                          self.volRatio, self.chargeRepeat, self.saltCon)
        self.pushoffResFileName = "res.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.pushoff".\
                                   format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                          self.volRatio, self.chargeRepeat, self.saltCon)

        self.pushoffLogFileName = "log.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.pushoff".\
                                  format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                         self.volRatio, self.chargeRepeat, self.saltCon)

        self.pushoffOutFileName = "out.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.pushoff".\
                                  format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                         self.volRatio, self.chargeRepeat, self.saltCon)

        if self.pushoff:

            try:
                self.in_pe_pushoff = open(self.in_pe_pushoffPath, 'w')
            except:
                pass

            content = "units lj\n" +\
                      "atom_style full\n" +\
                      "boundary p p p\n" +\
                      "processors * * *\n" +\
                      "read_data {0}\n".format(self.dataFileName) +\
                      "timestep 0.005\n" +\
                      "bond_style fene\n" +\
                      "bond_coeff 1 30.0 1.5 1.0 1.0\n" +\
                      "special_bonds fene\n" +\
                      "pair_style soft 1.0\n" +\
                      "pair_coeff * * 1.0 1.0\n" +\
                      "variable prefactor equal ramp(0,60)\n" +\
                      "fix 11 all adapt 1 pair soft a * * v_prefactor\n"  +\
                      "group polycation type 1 2\n" +\
                      "group polyanion type 3 4\n" +\
                      "group polyelectro type 1 2 3 4\n" +\
                      "group pcatcounter type 5\n" +\
                      "group pancounter type 6\n" +\
                      "group counterion type 5 6\n"
            
            if int(np.ceil(self.saltCon)) != 0:
                content += "group negsalt type 7\n" +\
                           "group possalt type 8\n" +\
                           "group salt type 7 8\n"
                    
            content += "velocity all create 1.0 864577 dist gaussian units box\n" +\
                       "variable r equal 'round(random(1, step+100, 48683))'\n" +\
                       "fix 1 all langevin 1.0 1.0 10.0 $r\n" +\
                       "fix 2 all nve\n" +\
                       "thermo 10000\n" +\
                       "thermo_style custom step atoms temp press etotal\n" +\
                       "thermo_modify lost warn flush yes\n" +\
                       "dump 1 all atom 10000 {0}\n".format(self.pushoffDumpFileName) +\
                       "restart 100000 {0}\n".format(self.pushoffResFileName) +\
                       "run 100000\n"

            self.in_pe_pushoff.write(content)
        
            self.in_pe_pushoff.close()


    def writeEquiInputFile(self):
        """
        write a equilibrium input file
        """
        self.in_pe_equiName = "la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.equi.{7}.txt".\
                              format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                     self.volRatio, self.chargeRepeat, self.saltCon, self.timestep)
        
        self.in_pe_equiPath = self.directory + self.in_pe_equiName

        self.equiDumpFileName = "dump.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.equi".\
                                format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                       self.volRatio, self.chargeRepeat, self.saltCon)

        self.equiResFileName = "res.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.equi".\
                               format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                      self.volRatio, self.chargeRepeat, self.saltCon)

        self.equiLogFileName = "log.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.equi.{7}.txt".\
                               format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                      self.volRatio, self.chargeRepeat, self.saltCon, self.timestep)

        self.equiOutFileName = "out.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.equi.{7}.txt".\
                               format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                      self.volRatio, self.chargeRepeat, self.saltCon, self.timestep)
        
        
        try:
            self.in_pe_equi = open(self.in_pe_equiPath, 'w')
        except:
            pass

        content =  "units lj\n" +\
                   "atom_style full\n" +\
                   "boundary p p p\n" +\
                   "processors * * *\n"

        if self.pushoff == False and self.first == True:
            if os.path.isfile(self.directory + self.pushoffResFileName + ".100000"):
                content += "read_restart " + self.pushoffResFileName + ".100000\n" +\
                           "reset_timestep 0\n"

            else:
                raise Exception("NO restart file founded!")

        if self.pushoff == False and self.first == False:
            if os.path.isfile(self.directory + self.equiResFileName + "." + str(self.pretimestep)):
                content += "read_restart " + self.equiResFileName + "." + str(self.pretimestep) + "\n"
                

            else:
                raise Exception("NO restart file founded!")


        content += "timestep 0.005\n" +\
                   "neigh_modify delay 4\n" +\
                   "pair_style none\n" +\
                   "pair_style lj/cut/coul/long/opt 1.122462048309373 10.0\n" +\
                   "pair_coeff 1*4 1*4 0.5 1.00 1.122462048309373\n"
                   
                   

        if int(np.ceil(self.saltCon)) != 0:
            content += "pair_coeff 1*4 5*8 0.5 0.75 0.841846536232029\n" +\
                       "pair_coeff 5*8 5*8 0.5 0.5 0.561231024154686\n"
        else:
            content += "pair_coeff 1*4 5*6 0.5 0.75 0.841846536232029\n" +\
                       "pair_coeff 5*6 5*6 0.5 0.50 0.561231024154686\n"
            
        content += "pair_modify shift yes\n" +\
                   "kspace_style pppm 1.0e-4\n" +\
                   "dielectric 0.576\n" +\
                   "bond_style fene\n" +\
                   "bond_coeff 1 30.0 1.5 1.0 1.0\n" +\
                   "special_bonds fene\n" +\
                   "group polycation type 1 2\n" +\
                   "group polyanion type 3 4\n" +\
                   "group polyelectro type 1 2 3 4\n" +\
                   "group pcatcounter type 5\n" +\
                   "group pancounter type 6\n" +\
                   "group counterion type 5 6\n"
        
        if int(np.ceil(self.saltCon)) != 0:
            content += "group negsalt type 7\n" +\
                       "group possalt type 8\n" +\
                       "group salt type 7 8\n"

        content += "variable r equal 'round(random(1, step+100, 48683))'\n" +\
                   "fix 1 all langevin 1.0 1.0 10.0 $r\n" +\
                   "fix 2 all nve\n" +\
                   "thermo 10000\n" +\
                   "thermo_style custom step atoms temp press etotal\n" +\
                   "thermo_modify lost warn flush yes\n" +\
                   "dump 1 all atom 10000 {0}\n".format(self.equiDumpFileName + "." + str(self.timestep)) +\
                   "restart 1000000 {0}\n".format(self.equiResFileName) +\
                   "run 2000000\n"

        self.in_pe_equi.write(content)
        self.in_pe_equi.close()


    def writePBSFile(self):
        """
        write a PBS file to be submitted on tbird.
        """
        if self.pushoff:
            self.pbsFileName = "mj.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.pushoff.pbs".\
                               format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                      self.volRatio, self.chargeRepeat, self.saltCon)
        else:
            self.pbsFileName = "mj.pe.la{0}.na{1}.lc{2}.nc{3}.rho{4}.r{5}.s{6}.equi.{7}.pbs".\
                               format(self.lenPa, self.numPa, self.lenPc, self.numPc, \
                                      self.volRatio, self.chargeRepeat, self.saltCon, self.timestep)

        self.pbsFilePath = self.directory + self.pbsFileName

        try:
            self.pbsFile = open(self.pbsFilePath, "w")
        except:
            pass

        content = "#!/bin/bash\n" +\
                  "#PBS -l walltime={0}:00:00\n".format(self.walltime) +\
                  "#PBS -l nodes={0}:ppn={1}\n".format(self.nodes, self.ppn) +\
                  "#PBS -j oe\n" +\
                  "cd $PBS_O_WORKDIR\n" +\
                  "date\n" +\
                  "mpirun -np {0} -machinefile $PBS_NODEFILE {1}".format(self.cores, self.lammps)

        if self.pushoff:
            content += " -in {0}".format(self.in_pe_pushoffName) +\
                       " -log {0}".format(self.pushoffLogFileName) +\
                       " > {0}\n".format(self.pushoffOutFileName)

        else:
            content += " -in {0}".format(self.in_pe_equiName) +\
                       " -log {0}".format(self.equiLogFileName) +\
                       " > {0}\n".format(self.equiOutFileName)


        self.pbsFile.write(content)
        self.pbsFile.close()

        if platform.node() == "tbird.phys.vt.edu":
            
            p = subprocess.Popen(['qsub', self.pbsFileName], cwd = self.directory)
            
                  

if __name__ == "__main__":
    job = Job(walltime = 240,
              nodes = 1,
              ppn = 16,
              lammps = "~/bin/lmp_tbird_15MAY2015",
              len_pan = 128,
              num_pan = 32,
              len_pcat = 32,
              num_pcat = 128,
              vol_ratio = 0.01,
              charge_fraction = 0.5,
              salt_concentration = 1.0,
              pretimestep = 0,
              timestep = 2000000,
              pushoff = False,
              first = False)
    
