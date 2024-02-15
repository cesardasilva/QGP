#!/usr/bin/env python3

from ROOT import TGraph, TCanvas, TF1
from array import array
from sys import argv

if len(argv) > 1:
    filename = argv[1]
else:
    filename = 'mydata.dat'

print('====== Reading', filename)
f = open(filename,'r')
lines = f.readlines()
x = array('f')
y = array('f')
for line in lines:
    if line[0][0] != '#':
        items = line.split()
        if items:
            x.append(float(items[0]))
            y.append(float(items[1]))
f.close()
print('======', len(x), 'data points read into arrays x and y')

print('====== Plotting data')
c1 = TCanvas('c1','Straight line fit example',10,10,600,400)
gr = TGraph(len(x),x,y)
gr.SetMarkerStyle(20)
gr.SetMarkerSize(1)
gr.Draw("AP")

print('====== Fitting data')
line = TF1('line', '[0]+[1]*x', min(x), max(x))
gr.Fit(line)
c1.Update()

input('====== Press Enter to exit...')
