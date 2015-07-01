#! /usr/bin/env python

import numpy as np

d='d.txt'
dat = np.loadtxt(d)
d2 = []
j=len(dat)
for i in range(1,j) :
    print dat[i]-dat[i-1]
