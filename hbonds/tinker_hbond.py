#! /usr/bin/env python

import h5py as h5
import sys
import numpy as np
import time
from optparse import OptionParser
import matplotlib.pyplot as plt
import pylab as P

def geometry( a1, a2, h, o ) :
    """
       
                 h-----o
                /
               /
       1------2
       
    """
    v21 = a1 - a2
    v21 /= np.linalg.norm(v21)
    v2h = h - a2
    v2h /= np.linalg.norm(v2h)
    vh2 = -v2h
    vho = o - h
    vho /= np.linalg.norm(vho)

    g12h = np.arccos(np.dot(v21, v2h))*180/np.pi
    g2ho = np.arccos(np.dot(vh2, vho))*180/np.pi
    return [ g12h, g2ho ]

parser = OptionParser()
parser.add_option("-F", "--h5file", dest="h5name", help=".h5 file")
parser.add_option("--atom1", dest="a1", help="Bond vector points AWAY from this atom")
parser.add_option("--atom2", dest="a2", help="Bond vector points TOWARDS this atom")
parser.add_option("--waterO", dest="ow", help="AMOEBA type for Water Oxygen.\tDefault: 247", default=247)
parser.add_option("--waterH", dest="oh", help="AMOEBA type for Water Hydrogen.\tDefault: 248", default=248)
parser.add_option("--r2H", dest="dist", help="Maximum distance between --atom2 and hydrogen.\tDefault: %.2f Angstroms"%(2.05+2*0.2), default=2.05+2*0.2)
parser.add_option("--a2HO", dest="a2ho", help="Minimum N--H-O angle.\tDefault: %.2f degrees"%(156-2*18), default=156-2*18)
parser.add_option("--a12H", dest="a12h", help="Minimum C-N--H angle.\tDefault: %.2f degrees"%(145-2*23), default=156-2*18)
parser.add_option("-o", "--out", dest="outpdf", help="Plot output name.\tDefault: out.pdf", default="out.pdf")
parser.add_option("-s", "--stop", dest="stop", help="Maximum number of frames to look at.\tDefault: None", default=None)

options, args = parser.parse_args()

a1 = int(options.a1) - 1
a2 = int(options.a2) - 1
ow = int(options.ow)
oh = int(options.oh)
dist = float(options.dist)
a2ho = float(options.a2ho)
a12h = float(options.a12h)
outpdf = options.outpdf
if outpdf[-4:] != ".pdf" : outpdf += ".pdf"

assert a1 >= 0
assert a2 >= 0
assert ow > 0
assert oh > 0
assert dist > 0
assert a2ho >= 0
assert a12h >= 0

f5 = h5.File(options.h5name, "r")
natoms = f5["t_atoms"].attrs.get("natoms")
nframes = f5["X"].attrs.get("nframes")
types = f5["t_atoms/type"][:]
# Get the index of water molecules
waters = []
water = []
for i in range(len(types)) :
    typei = types[i]
    if typei == ow :
        if len(water) == 3 :
            waters.append(water)
        water = [ i ]
    if typei == oh :
        water.append(i)

Rs = []
G12H = []
G2HO = []
W = []
fi = 0
for i in range(nframes) :
    if options.stop != None :
        if int(options.stop) == fi : break
    coords = f5["X/%i"%(i+1)]
    x1 = coords[a1]
    x2 = coords[a2]
    for each in waters :
        if np.linalg.norm(coords[each[0]] - x2) > dist + 2 : continue
        v1 = np.linalg.norm(coords[each[1]] - x2)
        v2 = np.linalg.norm(coords[each[2]] - x2)
        if v1 <= dist or v2 <= dist :
            if v1 < v2 :
                geo = geometry(x1, x2, coords[each[1]], coords[each[0]])
                if geo[0] >= a12h and geo[1] >= a2ho :
                    Rs.append(v1)
                    G12H.append(geo[0])
                    G2HO.append(geo[1])
                    W.append(each)
                    print "%i/%i"%(i,nframes), v1, geo[0], geo[1]
            else :
                geo = geometry(x1, x2, coords[each[2]], coords[each[0]])
                if geo[0] >= a12h and geo[1] >= a2ho :
                    Rs.append(v2)
                    G12H.append(geo[0])
                    G2HO.append(geo[1])
                    W.append(each)
                    print "%i/%i"%(i,nframes), v2, geo[0], geo[1]
    fi += 1


print len(Rs), len(Rs)/float(nframes)
outtxt = outpdf.replace(".pdf",".hbond")
ofile = open(outtxt, "w")
txtstring = "\t%-10i %-s\n"%(fi,options.h5name)
txtstring += ";%10s %20s %20s %20s\n"%("#HBonds", "<2H> Angstroms", "<12H> Degrees", "<2HO> Degrees")
txtstring += " %10i %10.2f%10.2f %10.2f%10.2f %10.2f%10.2f\n" %(len(Rs), np.average(Rs), np.std(Rs), np.average(G12H), np.std(G12H), np.average(G2HO), np.std(G2HO))
ofile.write(txtstring)
ofile.close()

def build_hist(ax, dat, title=None) :
    mu, sigma = np.average(dat), np.std(dat)
    n, bins, patches = ax.hist(dat, 50, normed=1, histtype='stepfilled')
    P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

    y = P.normpdf(bins, mu, sigma)
    l = ax.plot(bins, y, 'k--', linewidth=1.5)
    if title != None :
        ax.set_title(title)

fig = plt.figure()
ax1 = fig.add_subplot(223)
ax2 = fig.add_subplot(224, aspect='equal')
ax3 = fig.add_subplot(221)
ax4 = fig.add_subplot(222)
build_hist(ax1, Rs)
ax1.set_xlabel("N-H Distance ($\AA$)")
build_hist(ax3, G12H)
ax3.set_xlabel("C-N-H Angle (Degrees)")
build_hist(ax4, G2HO)
ax4.set_xlabel("N-H-O Angle (Degrees)")

ax2.hist2d(G12H, G2HO, bins=180/5, normed=True, range=[[0,180],[0,180]])
ax2.set_xlabel("C-N-H Angle")
ax2.set_ylabel("N-H-O Angle")
ax2.set_xticks(range(0,181,60))
ax2.set_yticks(range(0,181,60))
fig.subplots_adjust(hspace=0.35, wspace=0.3)
plt.savefig(outpdf, bbox_index='tight', pad_inches=0)

plt.close()
sys.exit()

fig = plt.figure()
ax1 = fig.add_subplot(223)
ax2 = fig.add_subplot(224, aspect='equal')
ax3 = fig.add_subplot(221)
ax4 = fig.add_subplot(222)
build_hist(ax1, Rs)
ax1.set_xlabel("N-H Distance ($\AA$)")
build_hist(ax3, G12H)
ax3.set_xlabel("C-N-H Angle (Degrees)")
build_hist(ax4, G2HO)
ax4.set_xlabel("N-H-O Angle (Degrees)")

ax2.hexbin(G12H, G2HO, extent=(0,180,0,180), xscale='linear', gridsize=180/5)
ax2.set_xlabel("C-N-H Angle")
ax2.set_ylabel("N-H-O Angle")
ax2.set_xticks(range(0,181,60))
ax2.set_yticks(range(0,181,60))
fig.subplots_adjust(hspace=0.35, wspace=0.3)
plt.savefig("hex+"+outpdf, bbox_index='tight', pad_inches=0)
