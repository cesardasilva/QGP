#!/usr/bin/env python3

from sys import argv
import numpy.random

npt = 15
yerr = 1.25
slope = 1.2
intc = -0.5

if len(argv) > 1:
    filename = argv[1]
    if len(argv) > 2:
        npt = int(argv[2])
else:
    filename = 'mydata.dat'
    npt = 15

f = open(filename,'w')

f.write('# Data for line fit test\n')
f.write('# intc = {}, slope = {}, yerr = {}\n'.format(intc,slope,yerr))

s = numpy.random.normal(0,yerr,npt)
for i in range(npt) :
    x = float(i+1)
    y = intc + slope*x + s[i]
    f.write('{}  {}\n'.format(x,y))

f.close()
print('====== Wrote', npt, 'points to', filename)
