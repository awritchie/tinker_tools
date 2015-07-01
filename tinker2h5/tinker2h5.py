#! /usr/bin/env python

import h5py as h5
import numpy as np
import sys
from optparse import OptionParser
import time

class Tinker2H5:
    def __init__ (self):
        self.ParseOptions()
        self.read_first_frame()
        self.read_input(self.xyzname,  "X")
        self.read_input(self.frcname,  "F")
        self.read_input(self.uindname, "U")
        self.read_input(self.velname,  "V")
        self.h5File.close()

    def read_first_frame(self) :
        if self.doAppend :
            self.natoms = self.h5File["t_atoms"].attrs.get("natoms")
            return
        print "\nReading first frame of %s"%self.xyzname
        natoms = 0
        index = 0
        names = []
        bonds = []
        types = []
        newFrame = False
        try :
            with open(self.xyzname,"rb") as file :
                for line in file :
                    l = line.split()
                    if newFrame :

                        index = 0

                    if len(l) == 1 :
                        if natoms != 0 : break
                        natoms = int(l[0])
                    else :
                        if self.isCoordLine(l) :
                            names.append( l[1] )
                            types.append( int(l[5]) )
                            bond = np.zeros(6) - 2
                            i = 0
                            for each in l[6:] :
                                bond[i] = int(each) - 1
                                i += 1
                            bonds.append(bond)
        except :
            sys.stderr.write("\nError reading %s, exiting.\n"%self.xyzname)
            sys.exit()
        names = np.array(names)
        bonds = np.array(bonds)
        types = np.array(types)
        self.natoms = len(names)
        grp = self.h5File.create_group("t_atoms")
        grp.create_dataset("name", data = names, dtype = "S10")
        grp.create_dataset("bond", data = bonds, maxshape = (self.natoms, 6))
        grp["bond"].attrs.create("NOTE", "Values of -2 are placeholders and do not indicate a molecular bond")
        grp.create_dataset("type", data = types)
        grp.attrs.create("natoms", data = self.natoms)
        print "Found %s atoms"%self.natoms
        return

    def read_input(self, inputName, type) :
        cstart = time.clock()
        wstart = time.time()
        if inputName == None : return
        print "\nReading %s"%inputName
        newFrame = False
        natoms = 0
        framei = 0
        if self.doAppend :
            framei = len(self.h5File[type])
            grp = self.h5File[type]
        else :
            grp = self.h5File.create_group(type)
            grp.attrs.create("nframes", data = 0)

        index = 0
        box = []
        X = []
        try :
            with open(inputName,"rb") as file :
                for line in file :
                    l = line.split()
                    if newFrame :
                        """ 
                                Make new array of coordinates
                                Start counting lines to get atom index
                                Index the frame number
                                Write the previous frame to file
                        """
                        if len(X) != 0 :
                            if self.natoms != index :
                                sys.stderr.write("Error! Frame %i has %i atoms, but we expected to have %i atoms.\n"%(framei, index, self.natoms))
                                return
                            grp.create_dataset("%i"%framei, data=X, compression=self.compress)
                        X = np.zeros(shape=(natoms,3))
                        index = 0
                        newFrame = False
                        framei += 1
                    if len(l) == 1 :
                        newFrame = True
                        natoms = int(l[0])
                    else :
                        if self.isCoordLine(l) :
                            X[index][0] = float(l[2].replace("D","E"))
                            X[index][1] = float(l[3].replace("D","E"))
                            X[index][2] = float(l[4].replace("D","E"))
                            index += 1
                        else :
                            framebox = np.array([ float(i) for i in l ])
                            box.append(framebox)

        except :
            sys.stderr.write("\nError reading %s, exiting.\n"%inputName)
            sys.exit()
        
        # Need to write out the last frame
        if len(X) != 0 :
            if self.natoms != index :
                sys.stderr.write("Error! Frame %i has %i atoms, but we expected to have %i atoms.\n"%(framei, index, self.natoms))
                return
            grp.create_dataset("%i"%framei, data=X, compression=self.compress)

        # Save box data per frame
        if len(box) > 0 :
            box = np.array(box)
            if self.doAppend :
                pBox = np.array(self.h5File["box"])
                tBox = np.concatenate((pBox, box))
                self.h5File["box"].resize((len(tBox),6))
                self.h5File["box"][:,:] = tBox
            else :
                self.h5File.create_dataset("box", data=box, maxshape = (None,6))
        grp.attrs.modify("nframes", framei)
        cend = time.clock()
        wend = time.time()
        if (wend - wstart) > 60*60*1.5 :
            print "Read %i frames (%.3f hr CPU time, %.3f hr Wallclock)"%(framei, (cend-cstart)/60/60, (wend-wstart)/60/60)
        elif (wend - wstart) > 60*1.5 :
            print "Read %i frames (%.3f min CPU time, %.3f min Wallclock)"%(framei, (cend-cstart)/60, (wend-wstart)/60/60)
        else :
            print "Read %i frames (%.3fs CPU time, %.3fs Wallclock)"%(framei, cend-cstart, wend-wstart)
    
    def isCoordLine(self,lineArray) :
        if len(lineArray) != 6 : return True
        if "." in lineArray[0] : return False
        if "." in lineArray[5] : return False
        return True
    
    def ParseOptions(self) :
        parser = OptionParser()
        parser.add_option("-F", "--h5file", dest="h5name", help="Output: .h5 file name", default="out.h5")
        parser.add_option("-x", "--xyz",  dest="xyzname", help="Input: coordinate {.xyz, .arc} file")
        parser.add_option("-f", "--frc",  dest="frcname", help="Input: force {.*f, .frc} file")
        parser.add_option("-u", "--uind", dest="uindname", help="Input: induced dipole {.*u, .uind} file")
        parser.add_option("-v", "--vel",  dest="velname", help="Input: velocity {.*v, .vel} file")
        parser.add_option("--append", dest="doAppend", help="Append to .h5 file.  Default: False", default=False, action="store_true")
        parser.add_option("--compress", dest="compress", help="Compress output .h5 file.  Expects either <gzip> or <lzf>.  Default None", default=None)
        options, args = parser.parse_args()

        self.h5name   = options.h5name
        self.xyzname  = options.xyzname
        self.frcname  = options.frcname
        self.uindname = options.uindname
        self.velname  = options.velname
        self.doAppend = options.doAppend
        self.compress = options.compress

        if self.doAppend :
            self.h5File = h5.File(self.h5name, "r+")
            self.compress = None
        else :
            if self.xyzname == None :
                sys.stderr.write("\nError: a coordinate file was not supplied.  Please supply a coordinate file with the --xyz option to create a .h5 file.\n")
                sys.exit()
            self.h5File = h5.File(self.h5name,"w")
            if self.compress != "gzip" and self.compress != "lzf" and self.compress != None :
                sys.stderr.write("\nWarning: <%s> unrecognized compression format.  Expected either <gzip> or <lzf>.  Defaulting to no compression.\n"%self.compress)
                self.compress = None




blah = Tinker2H5()
