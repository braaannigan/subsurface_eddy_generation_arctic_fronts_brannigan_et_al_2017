import os
import numpy as np

import genmit as gm
import eosmit
from arcf_mods import arcf_chdir as chd

def create_arcf_binary(z, dx, dy, nx, ny, f, H1, dh, dSh, Lf):

    nx = int(nx)
    ny = int(ny)
    dx = int(dx)
    dy = int(dy)
    f = float(f)
    H1 = int(H1)
    dh = float(dh)
    dSh = float(dSh)
    Lf = float(Lf)
    
    #Define the vertical grid
    z1=np.concatenate([-0.5*np.ones((50)) ,np.array([-0.6, -0.7, -0.8, -0.9]),
                       -1*np.ones((14)),-1.5*np.ones((32))] )
    x,y, bottom_depth = gm.gen_grid(dx,dy,nx,ny,z)
    gm.write_topography(dx,dy,nx,ny,bottom_depth)

    #Experiment parameters
    g = 9.81
    L = 0.5 * y[-1] #For mirrored domain
    Ys,Zs = np.meshgrid(y[0:int(0.5*ny)],z) #Short fields

    #Calculate the initial density field
    salt = 31. - 0.25* dSh * (1 + np.tanh( (Zs + H1) / dh)) * (1 + np.tanh( ( (-Ys + (2./3)*L)/Lf)-1  )  )
    salt = np.hstack( (np.fliplr(salt),salt) )
    temp = np.zeros( np.shape(salt) )

    Y,Z = np.meshgrid(y,z) #Short fields
    b = eosmit.buoy_calc(z,'lin',S = salt, Sref = 30)
    rho = 1025 -(1025/g)*b
    dz = np.abs(np.diff(z))
    dz = np.transpose(np.tile( (abs(np.diff(np.hstack( (0, z))))),(len(y),1)))

    eta = -(1./1025)*np.sum( (rho - 0.5 * (np.transpose(np.tile(rho[:,0] + rho[:,0.5*len(y)],(len(y),1))) )) * dz, axis = 0)


    u = np.zeros( np.shape(rho) )
    rhod = np.zeros( (len(z), len(y)+1))
    rhod[:,0:len(y)] = rho
    rhod[:,len(y)] = rho[:,-1]

    pre_factor = (g/f)*(1./1025)

    for i in np.arange(len(z)-1,-1,-1):
        u[i,:] = pre_factor * np.sum( dz[i:len(z),:] * np.diff(rhod[i:len(z),:])/dy,axis=0)

    temp[0,:] = 1e-4
    salt = np.transpose(np.tile(salt,(len(x),1,1)  ),(0,2,1))
    u = np.transpose(np.tile(u,(len(x),1,1)  ),(0,2,1))
    temp = np.transpose(np.tile(temp,(len(x),1,1)  ),(0,2,1))
    eta = np.tile(eta,(len(x),1) )
    #Add random perturbation
    amps=1e-3
    decay=H1
    SF = np.empty(np.shape(salt))
    for i in np.arange(0,len(z)):
        SF[:,:,i] = amps*np.random.rand(len(x),len(y))*np.exp(z[i]/decay)
    U = u + SF

    #Write the fields to binary
    temp = np.transpose( temp, (2,1,0))
    gm.write_temperature(dx,H1,temp)
    salt = np.transpose( salt, (2,1,0))
    gm.write_salt(dx,H1,salt)
    U = np.transpose( U, (2,1,0))
    gm.write_uvel(dx,H1,U)
    eta = np.transpose( eta, (1,0))
    gm.write_eta(dx,H1,eta)

