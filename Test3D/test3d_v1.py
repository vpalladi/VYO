#!/opt/local/bin/python2.7

import random
import math
import numpy as np
from array import array


from ROOT import gROOT, TFile, TH3D, TGraph

# parameters #
D = 0.5e-4    # um^2/s
L = 30.     # um
H = 30.     # um
W = 30.     # um
Nx = 10     # ! must be even !
Ny = 10     # ! must be even !
Nz = 10     # ! must be even !
dx = L/Nx   # um 
dy = H/Ny   # um
dz = W/Nz   # um
R =  15.    # um
dt = 3000.    # s
Nt = 1000  # numer of time steps

print "Simulating ",str(dt*Nt/60/60/24),"days"

C_contour = 0. # pg/um^3

# initial conditions
Ct0 = 0.021402   # pg/um^3
Ct0_sigma = Ct0*0.01

# Courant-Friedrichs-Lewy (CFL) restriction
Ratio_x = dt/(dx*dx)
Ratio_y = dt/(dy*dy)
Ratio_z = dt/(dz*dz)

if Ratio_x > (1/(2*D)) :
    print " !!! WARNING: Stability condiiton is not fulfilled !!!"
    print "dt/(dx*dx)",dt/(dx*dx)

    print "1/2*D",1/(2*D)
    print "2*D",2*D

print 'max value of dt with dx',str(dx),'is',str(dx*dx/(2*D))
print 'max value of dx with dt',str(dt),'is',str(math.sqrt(dt/(2*D)))

##############

## TGraph data ##

time = [0, 0.333, 2, 3, 5, 7, 10, 13, 16, 20, 21, 23, 27, 30, 34, 37, 42, 45]
time = array("d", [i*24. for i in time] )

concentration_perc = [0, 29, 37, 43, 48, 55, 63, 70, 74, 79, 80, 83, 82, 81, 82, 84, 79, 74]
concentration_perc = array("d", [(i/100.)+(i/100.)*(100-83)/100 for i in concentration_perc])

gData = TGraph(len(time), np.asarray(time), np.asarray(concentration_perc) )

#################

## numerical solution of dc/dt = D * V2c ##
def evolution(dt, dx, dy, dz, D, c_t_x_y_z, c_t_xp_y_z, c_t_xm_y_z, c_t_x_yp_z, c_t_x_ym_z, c_t_x_y_zp, c_t_x_y_zm) :
    A = D*dt/(dx*dx)
    B = D*dt/(dy*dy)
    C = D*dt/(dz*dz)

    Dx = c_t_xp_y_z + c_t_xm_y_z 
    Dy = c_t_x_yp_z + c_t_x_ym_z
    Dy = c_t_x_y_zp + c_t_x_y_zm
    
    c_tp_x_y_z = c_t_x_y_z * (1 - 2*A - 2*B - 2*C) + Dx * A + Dy * B + Dy * C

    return c_tp_x_y_z
############################################

# out file 
FileOutRoot = TFile('Out.root', 'recreate')

# position and concentration matrixes 
P = [[[[0.,0., 0.] for iz in range(Nz+2)] for iy in range(Ny+2)] for ix in range(Nx+2)]
C = [[[C_contour for iz in range(Nz+2)] for iy in range(Ny+2)] for ix in range(Nx+2)]

# graphs
gMassEvolution = TGraph()
gMassRelease = TGraph()

# generate the map at t0
Cmap_t0 = TH3D("Cmap_t0", "Cmap_t0",
               Nx+2, (-Nx/2-1)*dx+dx/2, (Nx/2+1)*dx+dx/2,
               Ny+2, (-Ny/2-1)*dy+dy/2, (Ny/2+1)*dy+dy/2,
               Nz+2, (-Nz/2-1)*dz+dz/2, (Nz/2+1)*dz+dz/2 )

Total_amount_solute = 0 # g

print " > Generating the map at t0"
for ix in range(0, Nx+2) :
    for iy in range(0, Ny+2) :
        for iz in range(0, Nz+2) :
        
            X = (ix-Nx/2-1)*dx
            Y = (iy-Ny/2-1)*dy
            Z = (iz-Nz/2-1)*dz
            r = math.sqrt( X*X + Y*Y + Z*Z )
        
            P[ix][iy][iz] = [X, Y, Z]
            
            if r<R :
                value = random.gauss(Ct0, Ct0_sigma)
                if value<0 :
                    value = 0 
                C[ix][iy][iz] = value
                
            Px = P[ix][iy][iz][0]
            Py = P[ix][iy][iz][1]
            Pz = P[ix][iy][iz][2]

            iBinX = Cmap_t0.GetXaxis().FindBin( Px )
            iBinY = Cmap_t0.GetYaxis().FindBin( Py )
            iBinZ = Cmap_t0.GetZaxis().FindBin( Pz )
            Cmap_t0.SetBinContent( iBinX, iBinY, iBinZ, C[ix][iy][iz] )

            Total_amount_solute += dx*dy*dz*C[ix][iy][iz]

gMassEvolution.SetPoint( gMassEvolution.GetN(), 0, Total_amount_solute )
print ' > initial amount of solute:',Total_amount_solute,'pg'  
Cmap_t0.Write("T0_C")

## time evolution ##
print " > Evolving the matrix"
for t in range(1, Nt) :

    # counter
    if int(t%(Nt*0.1)) == 0 :
        print str(int(10*t/(Nt*0.1)))+"% completed"

    C_tmp = [[[C_contour for iz in range(Nz+2)] for iy in range(Ny+2)] for ix in range(Nx+2)]

    # LOOP over all the pixels
    for ix in range(1, Nx+1) :
        for iy in range(1, Ny+1) :
            for iz in range(1, Nz+1) :

                r = math.sqrt(P[ix][iy][iz][0]*P[ix][iy][iz][0] + P[ix][iy][iz][1]*P[ix][iy][iz][1])
                if r<R : 
                    C_tmp[ix][iy][iz] =  evolution(dt, dx, dy, dz, D,
                                                   C[ix][iy][iz],
                                                   C[ix+1][iy][iz], C[ix-1][iy][iz],
                                                   C[ix][iy+1][iz], C[ix][iy-1][iz],
                                                   C[ix][iy][iz+1], C[ix][iy][iz-1] )

    Cmap = TH3D("Cmap_t"+str(t), "Cmap_t"+str(t),
                Nx+2, (-Nx/2-1)*dx+dx/2, (Nx/2+1)*dx+dx/2,
                Ny+2, (-Ny/2-1)*dy+dy/2, (Ny/2+1)*dy+dy/2,
                Nz+2, (-Nz/2-1)*dz+dz/2, (Nz/2+1)*dz+dz/2 )

    IntegralTm = 0
    IntegralT  = 0
    
    for ix in range(0, Nx+2) :
        for iy in range(0, Ny+2) :
            for iz in range(0, Nz+2) :

                IntegralTm += dx*dy*dz*C[ix][iy][iz] 
                
                C[ix][iy][iz] = C_tmp[ix][iy][iz]

                IntegralT += dx*dy*dz*C[ix][iy][iz] 

                Px = P[ix][iy][iz][0]
                Py = P[ix][iy][iz][1]
                Pz = P[ix][iy][iz][2]

                iBinX = Cmap_t0.GetXaxis().FindBin( Px )
                iBinY = Cmap_t0.GetYaxis().FindBin( Py )
                iBinZ = Cmap_t0.GetZaxis().FindBin( Pz )

                Cmap.SetBinContent( iBinX, iBinY, iBinZ, C[ix][iy][iz] )

    gMassEvolution.SetPoint( gMassEvolution.GetN(), t*dt/3600, IntegralT )
    gMassRelease.SetPoint( gMassRelease.GetN(), t*dt/3600, (Total_amount_solute-IntegralT)/Total_amount_solute )

#    Cmap.Write( 'Ct_'+str(t*dt) )

print "Completed!"

gMassEvolution.Write("Mass_Evolution")
gMassRelease.Write("Mass_Release")
gData.Write("YData")

FileOutRoot.Close()
