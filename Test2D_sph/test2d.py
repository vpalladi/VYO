#!/opt/local/bin/python2.7

import random
import math
import numpy as np

from ROOT import gROOT, TFile, TH2D, TGraph

# parameters #
D = 0.5e-4   # cm^2/s
R = 30.      # um
Nr = 10      # ! must be even !
Nteta = 10   # ! must be even !
dr = R/Nr    # um 
dteta = math.pi/Nteta # um

dt = 600.    # s
Nt = 10000   # numer of time steps

print 'D     = ',D    ,' cm^2/s'
print 'R     = ',R    ,' um    '
print 'Nr    = ',Nr            
print 'Nteta = ',Nteta         
print 'dr    = ',dr   ,' um    '
print 'dteta = ',dteta,' um    '
print 'dt    = ',dt   ,' s     '
print 'Nt    = ',Nt   

print "Simulating ",str(dt*Nt/60),"minutes"

C_contour = 0. # g/cm^3

Ct0 = 0.021402   # pg/um^3
Ct0_sigma = Ct0*0.01

### # Courant-Friedrichs-Lewy (CFL) restriction
### Ratio = dt/(dx*dx)
### 
### if Ratio > (1/(2*D)) :
###     print " !!! WARNING: Stability condiiton is not fulfilled !!!"
###     print "dt/(dx*dx)",dt/(dx*dx)
###     print "1/2*D",1/2*D

##############

## numerical solution of dc/dt = D * V2c ##
def evolution( dt, dr, r, dteta, D, c_t_r_teta, c_t_rp_y, c_t_rm_teta, c_t_r_tetap, c_t_r_tetam) :
    a = 1 - D * dt * ( 1/(r*dr) + 2/(dr*dr) + 2/(dteta*dteta) )
    b = 1/(r*dr)
    c = 1/(dr*dr)
    d = 1/(r*r*dteta*dteta)
    Dr    = c_t_rp_teta + c_t_rm_teta
    Dteta = c_t_r_tetap + c_t_r_tetam
    
    c_tp_r_teta = c_t_r_teta * a + c_t_rp_teta * b + Dr * c + Dteta * d  
 
    return c_tp_r_teta
#########################################

# out file 
FileOutRoot = TFile('Out.root', 'recreate')

# position and concentration matrixes 
P = [[[0.,0.] for iteta in range(Nteta)] for ir in range(Nr+2)]
C = [[C_contour for iteta in range(Nteta)] for ir in range(Nr+2)]

# generate the map at t0
Cmap_t0 = TH2D("Cmap_t0", "Cmap_t0",
               Nteta, 0, Nteta*dteta,
               Nr+2, 0, (Nr+2)*dr )

print " > Generating the map at t0"
for ir in range(0, Nr+2) :
    for iteta in range(0, Nteta) :
    
        teta = iteta * dteta
        r = ir * dr
        print ir, iteta
        print r, teta
        
        P[ir][iteta] = [r, teta]

        if r<=R :
            value = random.gauss(Ct0, Ct0_sigma)
            if value<0 :
                value = 0 
            C[ir][iteta] = value
                
        Pr = P[ir][iteta][0]
        Pteta = P[ir][iteta][1]

        iBinX = Cmap_t0.GetXaxis().FindBin( Pteta )
        iBinY = Cmap_t0.GetYaxis().FindBin( Pr )
        Cmap_t0.SetBinContent( iBinX, iBinY, C[ir][iteta] )
                
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
