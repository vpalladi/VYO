#!/opt/local/bin/python2.7

import random
import math
import numpy as np

from ROOT import gROOT, TFile, TH2D, TGraph

# parameters #
D = 8e-8   # cm^2/s
L = 1.     # cm
H = 1.     # cm
Nx = 10    # ! must be even !
Ny = 10    # ! must be even !
dx = L/Nx   # cm 
dy = H/Ny   # cm
R = 0.4      # cm
dt = 600.  # s
Nt = 10000   # numer of time steps

print "Simulating ",str(dt*Nt/60),"minutes"

C_contour = 0. # g/cm^3

Ct0 = 0.1   # g/cm^3
Ct0_sigma = Ct0*0.01

# Courant-Friedrichs-Lewy (CFL) restriction
Ratio = dt/(dx*dx)

if Ratio > (1/(2*D)) :
    print " !!! WARNING: Stability condiiton is not fulfilled !!!"
    print "dt/(dx*dx)",dt/(dx*dx)
    print "1/2*D",1/2*D

##############

## numerical solution of dc/dt = D * V2c ##
def evolution(dt, dx, dy, D, c_t_x_y, c_t_xp_y, c_t_xm_y, c_t_x_yp, c_t_x_ym) :
    A = D*dt/(dx*dx)
    B = D*dt/(dy*dy)
    # C = D*dt/(dz*dz)
    Dx = c_t_xp_y + c_t_xm_y 
    Dy = c_t_x_yp + c_t_x_ym
    # Dy = c_t_x_y_zp + c_t_x_y_zm 
    c_tp_x_y = c_t_x_y * (1 - 2*A - 2*B) + Dx * A + Dy * B

    return c_tp_x_y
############################################

# out file 
FileOutRoot = TFile('Out.root', 'recreate')

# position and concentration matrixes 
P = [[[0.,0.] for iy in range(Ny+2)] for ix in range(Nx+2)]
C = [[C_contour for iy in range(Ny+2)] for ix in range(Nx+2)]

# generate the map at t0
Cmap_t0 = TH2D("Cmap_t0", "Cmap_t0",
               Nx+2, (-Nx/2-1)*dx+dx/2, (Nx/2+1)*dx+dx/2,
               Ny+2, (-Ny/2-1)*dy+dy/2, (Ny/2+1)*dy+dy/2 )


print " > Generating the map at t0"
for ix in range(0, Nx+2) :
    for iy in range(0, Ny+2) :

        X = (ix-Nx/2-1)*dx
        Y = (iy-Ny/2-1)*dy
        r = math.sqrt( X*X + Y*Y )
        
        P[ix][iy] = [X, Y]

        if r<R :
            value = random.gauss(Ct0, Ct0_sigma)
            if value<0 :
                value = 0 
            C[ix][iy] = value
                
        Px = P[ix][iy][0]
        Py = P[ix][iy][1]

        iBinX = Cmap_t0.GetXaxis().FindBin( Px )
        iBinY = Cmap_t0.GetYaxis().FindBin( Py )
        Cmap_t0.SetBinContent( iBinX, iBinY, C[ix][iy] )
                
Cmap_t0.Write("T0_C")

## time evolution ##

print " > Evolving the matrix"
for t in range(1, Nt) :
    if int(t%(Nt*0.1)) == 0 :
        print str(int(10*t/(Nt*0.1)))+"% completed"
        
    C_tmp = [[C_contour for iy in range(Ny+2)] for ix in range(Nx+2)]

    for ix in range(1, Nx+1) :
        for iy in range(1, Ny+1) :
            r = math.sqrt(P[ix][iy][0]*P[ix][iy][0] + P[ix][iy][1]*P[ix][iy][1])
            if r<R : 
                C_tmp[ix][iy] =  evolution(dt, dx, dy, D,
                                           C[ix][iy],
                                           C[ix+1][iy], C[ix-1][iy],
                                           C[ix][iy+1], C[ix][iy-1]) 

    Cmap = TH2D("Cmap_t"+str(t), "Cmap_t"+str(t),
                Nx+2, (-Nx/2-1)*dx+dx/2, (Nx/2+1)*dx+dx/2,
                Ny+2, (-Ny/2-1)*dy+dy/2, (Ny/2+1)*dy+dy/2 )
    
    for ix in range(0, Nx+2) :
        for iy in range(0, Ny+2) :

            C[ix][iy] = C_tmp[ix][iy]

            Px = P[ix][iy][0]
            Py = P[ix][iy][1]

            iBinX = Cmap_t0.GetXaxis().FindBin( Px )
            iBinY = Cmap_t0.GetYaxis().FindBin( Py )

            Cmap.SetBinContent( iBinX, iBinY, C[ix][iy] )
            
    Cmap.Write( 'Ct_'+str(t*dt) )
print "Completed!"
    
FileOutRoot.Close()
