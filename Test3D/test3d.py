#!/opt/local/bin/python2.7

import random
import math
import numpy as np
from array import array

from ROOT import gROOT, TFile, TH1D, TH3D, TGraph, TGraphErrors


# parameters #
D = 1e-4      # um^2/s
D_crit = 0.
L = 28.       # um
H = 28.       # um
W = 28.       # um
Nx = 20       # ! must be even ! # in the majority of the following the index goes from 0 to Nx+2 in order to consider a pixel as the external medium 
Ny = 20       # ! must be even ! # in the majority of the following the index goes from 0 to Nx+2 in order to consider a pixel as the external medium 
Nz = 20       # ! must be even ! # in the majority of the following the index goes from 0 to Nx+2 in order to consider a pixel as the external medium 
dx = L/Nx     # um 
dy = H/Ny     # um
dz = W/Nz     # um
R =  14.0      # um
pixel_meanlife = 144. # hours
# polimer density 1.4268 +- 0.0003 g/cm^3

if R > L :
    print " !!! WARNING sphere radius must be bigger than volume !!!"
    exit
    
dt = 3600.     # s
Nt = 72      # numer of time steps



print "Simulating ",str(dt*Nt/60/60/24),"days"

C_contour = 0. # pg/um^3

## initial conditions

# concentration 15ug/mg (before 25/4/2015) >>>> Ct0 = 0.021402 pg/um^3
# concentration 10ug/mg (25/4/2015) >>>> Ct0 = 0.021402*0.66 pg/um^3
Ct0 = 0.021402*0.66 ## pg/um^3 
Ct0_sigma = Ct0*0.01

## Courant-Friedrichs-Lewy (CFL) restriction
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


## numerical solution of dc/dt = D * V2c ##

def evolution(dt, dx, dy, dz, D, D_crit, nx_no_deg, ny_no_deg, nz_no_deg, nx_tot, ny_tot, nz_tot, c_t_x_y_z, c_t_xp_y_z, c_t_xm_y_z, c_t_x_yp_z, c_t_x_ym_z, c_t_x_y_zp, c_t_x_y_zm) :
    D_eff_x = D + D_crit * (1-nx_no_deg/nx_tot)
    D_eff_y = D + D_crit * (1-ny_no_deg/ny_tot)
    D_eff_z = D + D_crit * (1-nz_no_deg/nz_tot)

    A = D_eff_x*dt/(dx*dx)
    B = D_eff_y*dt/(dy*dy)
    C = D_eff_z*dt/(dz*dz)

    Dx = float(c_t_xp_y_z) + float(c_t_xm_y_z) 
    Dy = float(c_t_x_yp_z) + float(c_t_x_ym_z)
    Dy = float(c_t_x_y_zp) + float(c_t_x_y_zm)

    c_tp_x_y_z = float(c_t_x_y_z) * (1 - 2*A - 2*B - 2*C) + Dx * A + Dy * B + Dy * C

    return c_tp_x_y_z

############################################

## class pixel ##

class Cpixel :
    
    def __init__(self, pixel_meanlife_h) :
        self.__mean_life = pixel_meanlife_h
        self.time_of_degradation_h = random.expovariate(1./self.__mean_life)
        print  self.time_of_degradation_h
    dx = 0.
    dy = 0.
    dz = 0.
    center = [0., 0., 0.] # [x, y, z]
    r = math.sqrt( center[0]*center[0] + center[1]*center[1] + center[2]*center[2] )
    C_0 = 0.
    C = 0.
    D = 0.
#    is_degradated = 0

################# evolve cannot work as it is: the evolution of the matrix is such that once one pixel i evolved the neighborhood pixels will use the value at t and not a t-1
#    def evolve(self, dt, c_t_xp_y_z, c_t_xm_y_z, c_t_x_yp_z, c_t_x_ym_z, c_t_x_y_zp, c_t_x_y_zm) :
#        if self.is_degradated == 0 :
#            if self.__first_iteration == 1 :
#                self.__first_iteration = 0
#                self.C = self.C_0
#                self.C = evolution(dt, self.dx, self.dy, self.dz, self.D, self.C_0,
#                                   c_t_xp_y_z, c_t_xm_y_z,
#                                   c_t_x_yp_z, c_t_x_ym_z,
#                                   c_t_x_y_zp, c_t_x_y_zm )
#            else :
#                self.C = evolution(dt, self.dx, self.dy, self.dz, self.D, self.C,
#                                   c_t_xp_y_z, c_t_xm_y_z,
#                                   c_t_x_yp_z, c_t_x_ym_z,
#                                   c_t_x_y_zp, c_t_x_y_zm )
                
    def is_degradate(self, T_h) :
        if T_h >= self.time_of_degradation_h :
            return 1
        else :
            return 0

    __meanlife = 0.
    __first_iteration = 1
    time_of_degradation_h = 0. 
############################################


## TGraph data ##

time = [0, 0.333, 2, 3, 5, 7, 10, 13, 16, 20, 21, 23, 27, 30, 34, 37, 42, 45]
time = array( "d", [i*24. for i in time] )

concentration_perc = [0, 29, 37, 43, 48, 55, 63, 70, 74, 79, 80, 83, 82, 81, 82, 84, 79, 74]
concentration_perc = array("d", [(i/100.)+(i/100.)*(100-83)/100 for i in concentration_perc])

gData = TGraph(len(time), np.asarray(time), np.asarray(concentration_perc) )


time_2 = [1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 36, 38, 40, 43, 44.5, 46, 48] # hours
time_err_2 = array( "d", [0.0833 for i in time_2] )

concentration_perc_2 = [4.75, 6.53, 7.94, 8.87, 11.11, 13.02, 13.92, 16.91, 17.94, 18.82, 19.72, 20.42, 20.84, 21.15, 21.90, 22.32, 24.82, 25.33, 25.84, 26.50, 26.82, 27.26, 27.57]
concentration_perc_2 = array( "d", [(i/100.) for i in concentration_perc_2])

concentration_err_2 = [0.8, 0.8, 0.9, 1.0, 1.1, 1.2, 1.2, 0.3, 0.3, 0.4, 0.4, 0.2, 0.4, 0.3, 0.5, 0.4, 0.5, 0.4, 0.5, 0.4, 0.5, 0.4, 0.3]
concentration_err_2 = array( "d", [(i/100.) for i in concentration_err_2])

gData_2 = TGraphErrors(len(time_2), np.asarray(time_2), np.asarray(concentration_perc_2), np.asarray(time_err_2), np.asarray(concentration_err_2) )

###############################################

# out file

FileOutName = 'Out_R'+str(R)+'.root'
FileOutRoot = TFile(FileOutName, 'recreate')

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

# define the matrix of pixels
MATRIX = [[[Cpixel(pixel_meanlife) for i in range(Nx+2)] for i in range(Ny+2)] for i in range(Nz+2)]

hDegradation_time_distribution = TH1D('degradation_time_distribution', 'degradation_time_distribution',
                                      100, 0, pixel_meanlife*4)

print " > Generating the map at t0"
for ix in range(0, Nx+2) :
    for iy in range(0, Ny+2) :
        for iz in range(0, Nz+2) :
            hDegradation_time_distribution.Fill( MATRIX[ix][iy][iz].time_of_degradation_h )
            X = (ix-Nx/2-1)*dx
            Y = (iy-Ny/2-1)*dy
            Z = (iz-Nz/2-1)*dz
            r = math.sqrt( X*X + Y*Y + Z*Z )
        
            MATRIX[ix][iy][iz].center = [X, Y, Z]
            MATRIX[ix][iy][iz].D = D
            MATRIX[ix][iy][iz].dx = dx
            MATRIX[ix][iy][iz].dy = dy
            MATRIX[ix][iy][iz].dz = dz

            value = 0.
            if r<R :
                value = random.gauss(Ct0, Ct0_sigma)
                if value<0 :
                    value = 0
   
            MATRIX[ix][iy][iz].C_0 = value
            MATRIX[ix][iy][iz].C = value

            
            iBinX = Cmap_t0.GetXaxis().FindBin( MATRIX[ix][iy][iz].center[0] )
            iBinY = Cmap_t0.GetYaxis().FindBin( MATRIX[ix][iy][iz].center[1] )
            iBinZ = Cmap_t0.GetZaxis().FindBin( MATRIX[ix][iy][iz].center[2] )
            Cmap_t0.SetBinContent( iBinX, iBinY, iBinZ, MATRIX[ix][iy][iz].C_0 )

            Total_amount_solute += dx*dy*dz*MATRIX[ix][iy][iz].C_0

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

    for ix in range(1, Nx+1) :
        for iy in range(1, Ny+1) :
            for iz in range(1, Nz+1) :
                C_tmp[ix][iy][iz] = MATRIX[ix][iy][iz].C

    # LOOP over all the pixels
    for ix in range(1, Nx+1) :
        for iy in range(1, Ny+1) :
            for iz in range(1, Nz+1) :
                
                r = MATRIX[ix][iy][iz].r
                
                if r<R :
                    MATRIX[ix][iy][iz].C = evolution(dt, dx, dy, dz, D, D_crit,
                                                     1, 1, 1,
                                                     1, 1, 1,
                                                     C_tmp[ix][iy][iz],
                                                     C_tmp[ix+1][iy][iz], C_tmp[ix-1][iy][iz],
                                                     C_tmp[ix][iy+1][iz], C_tmp[ix][iy-1][iz],
                                                     C_tmp[ix][iy][iz+1], C_tmp[ix][iy][iz-1] )

                #MATRIX[ix][iy][iz].degradate(t*dt/3600.)
                    
    Cmap = TH3D("Cmap_t"+str(t), "Cmap_t"+str(t),
                Nx+2, (-Nx/2-1)*dx+dx/2, (Nx/2+1)*dx+dx/2,
                Ny+2, (-Ny/2-1)*dy+dy/2, (Ny/2+1)*dy+dy/2,
                Nz+2, (-Nz/2-1)*dz+dz/2, (Nz/2+1)*dz+dz/2 )
    IntegrityMap = TH3D("IntegrityMap"+str(t), "IntegrityMap"+str(t),
                        Nx+2, (-Nx/2-1)*dx+dx/2, (Nx/2+1)*dx+dx/2,
                        Ny+2, (-Ny/2-1)*dy+dy/2, (Ny/2+1)*dy+dy/2,
                        Nz+2, (-Nz/2-1)*dz+dz/2, (Nz/2+1)*dz+dz/2 )

    IntegralTm = 0
    IntegralT  = 0
    
    for ix in range(0, Nx+2) :
        for iy in range(0, Ny+2) :
            for iz in range(0, Nz+2) :

                IntegralTm += dx*dy*dz*C_tmp[ix][iy][iz] 
                IntegralT  += dx*dy*dz*MATRIX[ix][iy][iz].C
                
                iBinX = Cmap_t0.GetXaxis().FindBin( MATRIX[ix][iy][iz].center[0] )
                iBinY = Cmap_t0.GetYaxis().FindBin( MATRIX[ix][iy][iz].center[1] )
                iBinZ = Cmap_t0.GetZaxis().FindBin( MATRIX[ix][iy][iz].center[2] )

                Cmap.SetBinContent( iBinX, iBinY, iBinZ, C[ix][iy][iz] )
                if MATRIX[ix][iy][iz].is_degradate(t*dt/3600) :
                    IntegrityMap.SetBinContent( iBinX, iBinY, iBinZ, 0 )
                else :
                    IntegrityMap.SetBinContent( iBinX, iBinY, iBinZ, 1 )
    gMassEvolution.SetPoint( gMassEvolution.GetN(), t*dt/3600, IntegralT )
    gMassRelease.SetPoint( gMassRelease.GetN(), t*dt/3600, (Total_amount_solute-IntegralT)/Total_amount_solute )

    IntegrityMap.Write( 'IntegrityMap_t'+str(t*dt) )
#    Cmap.Write( 'Ct_'+str(t*dt) )

print "Completed!"

gMassEvolution.Write("Mass_Evolution")
gMassRelease.Write("Mass_Release")
gData.Write("YData")
gData_2.Write("YData_april2015")
hDegradation_time_distribution.Write()


FileOutRoot.Close()
