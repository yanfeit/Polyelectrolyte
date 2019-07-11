#!/bin/usr/python
#

docstr = """
Split the large dump file into small one.
usage: SplitDump.py sys.argv[1] sys.argv[2]
"""

import sys
import numpy as np
import os

class SplitDump():

    def __init__(self, iname, oname):

        f = open(iname, "r")

        while True:
            line = f.readline()
            if line == '':
                break

            elif line == "ITEM: TIMESTEP\n":
                ofile = open(oname, "w")
                ofile.write(line)  # ITEM
                line = f.readline() # 0
                time = int(line)
                ofile.write(line) # 0
                line = f.readline() # ITEM: NUMBER OF ATOMS
                ofile.write(line) #
                line = f.readline() # 6600
                natoms = int(line)
                ofile.write(line)
                for i in range(natoms + 5):
                    line = f.readline()
                    ofile.write(line)
                    ofile.close()
                    os.rename(oname, sys.argv[2]+"{0}".format(time))

        
        f.close()


if __name__ == "__main__":
    sample = SplitDump(sys.argv[1], sys.argv[2])


