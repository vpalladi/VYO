#!/opt/local/bin/python2.7

import random
from ROOT import gROOT, TFile, TH2Poly
import numpy as np

D = 8e-8   # cm^2/s
L = 1.     # cm
H = 1.     # cm
Nx = 50
Ny = 50   
dx = L/Nx   # cm 
dy = H/Ny   # cm
dt = 600.  # s
Nt = 1000   # numer of time steps

C_contour = 0. # g/cm^3

Ct0 = 0.1   # g/cm^3
Ct0_sigma = Ct0*0.01

def evolution(dt, dx, dy, D, c_t_x_y, c_t_xp_y, c_t_xm_y, c_t_x_yp, c_t_x_ym) :
    A = D*dt/(dx*dx)
    B = D*dt/(dy*dy)
    Dx = c_t_xp_y + c_t_xm_y 
    Dy = c_t_x_yp + c_t_x_ym 
    c_tp_x_y = c_t_x_y * (1 - 2*A - 2*B) + Dx * A + Dy * B

    return c_tp_x_y

FileOutRoot = TFile('Out.root', 'recreate')

# position
P = [[[0.,0.] for iy in range(Ny+2)] for ix in range(Nx+2)]
# concentration
C = [[C_contour for iy in range(Ny+2)] for ix in range(Nx+2)]

Cmap_t0 = TH2Poly()

for ix in range(0, Nx+2) :
    for iy in range(0, Ny+2) :

        P[ix][iy] = [(ix-1)*dx, (iy-1)*dy]
        
        if ix!=0 and iy!=0 and ix!=Nx+1 and iy!=Ny+1 :
            value = random.gauss(Ct0, Ct0_sigma)
            if value<0 :
                value = 0
            C[ix][iy] = value
                    
        Px = P[ix][iy][0]
        Py = P[ix][iy][1]
      
        Cmap_t0.AddBin(Px-dx/2, Py-dy/2, Px+dx/2, Py+dy/2)
        Cmap_t0.SetBinContent( Cmap_t0.GetNumberOfBins(), C[ix][iy] )

Cmap_t0.Write("T0_C")

## time evolution ##
for t in range(1, Nt) :
    C_tmp = [[C_contour for iy in range(Ny+2)] for ix in range(Nx+2)]

    for ix in range(1, Nx+1) :
        for iy in range(1, Ny+1) :
            C_tmp[ix][iy] =  evolution(dt, dx, dy, D,
                                       C[ix][iy],
                                       C[ix+1][iy], C[ix-1][iy],
                                       C[ix][iy+1], C[ix][iy-1]) 

    Cmap = TH2Poly()
    
    for ix in range(0, Nx+2) :
        for iy in range(0, Ny+2) :

            C[ix][iy] = C_tmp[ix][iy]

            Px = P[ix][iy][0]
            Py = P[ix][iy][1]
            Cmap.AddBin(Px-dx/2, Py-dy/2, Px+dx/2, Py+dy/2)
            Cmap.SetBinContent( Cmap.GetNumberOfBins(), C[ix][iy] )
            
    Cmap.Write( 'Ct_'+str(t*dt) )

    
FileOutRoot.Close()
