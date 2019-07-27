#!/opt/local/bin/python2.7

import random
from ROOT import gROOT, TCanvas, TGraph, TFile
import numpy as np

D = 8e-8 # cm^2/s
L = 1.     # cm
N = 100   # 
dx = L/N  # cm 
dt = 500.    # s

C_0 = 0. # g/cm^3
C_L = 0.   # g/cm^3

Ct0 = 0.1   # g/cm^3

def evolution(dt, dx, D, c_i_j, c_i_jp, c_i_jm) :
    A = dt/(dx*dx)    
    B = c_i_jp + c_i_jm
    C = 1 - 2 * D * A
    c_ip_j = D * A * B + c_i_j * C

    return c_ip_j

print 'L  ', L
print 'D  ', D
print 'N  ', N
print 'dt ', dt
print 'dx ', dx

print str(dt/(dx*dx))+'<'+str(1/D)

C = [Ct0]
X = [0]

for i in range(1, N-1) :
    C.append( random.gauss(Ct0, Ct0*0.001) )
    X.append( X[-1]+dx )
    
C.append( C_L ) 
X.append( L )

FileOutRoot = TFile('Out.root', 'recreate')

for t in range(0, 10000) :
    C_tp = [C_0]
    for i in range(1, N-1) :
        C_tp.append( evolution( dt, dx, D, C[i], C[i+1], C[i-1] ) )
    C_tp.append( C_L )
    
    g = TGraph( N, np.asarray(X), np.asarray(C_tp) )
    g.Write( 'Ct_'+str(t) )
    
    for i in range(0, N) :
        C[i] = C_tp[i]
        
FileOutRoot.Close()
