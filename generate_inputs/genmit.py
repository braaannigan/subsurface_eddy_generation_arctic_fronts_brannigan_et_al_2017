from __future__ import print_function

import numpy as np

#Write binary files with big endian entries
def writefield(file_name,data):
    import sys
    print(('write to file: '+ file_name))
    print((np.shape(data)))
    if sys.byteorder == 'little': data.byteswap(True)
    fid = open(file_name,"wb")
    data.tofile(fid)
    fid.close()

def gen_grid(dx,dy,nx,ny,z):
    x=np.arange(0.5*dx,nx*dx+0.5*dx,dx)
    y=np.arange(0.5*dy,ny*dy+0.5*dy,dy)
    bottom_depth = z[-1]*np.ones( (nx,ny) )
    return x, y, bottom_depth

def write_topography(dx,dy,nx,ny,bottom_depth):
    file_name = 'dx' + str(dx) + 'nx' + str(nx) + 'ny' + str(ny) +'topog.bin'
    writefield(file_name, bottom_depth)

def write_temperature(dx,H1,temp):
    file_name = 'dx' + str(dx) + 'H1' +str(H1) + 'PoTemp.bin'
    writefield(file_name, temp)

def write_salt(dx,H1,salt):
    file_name = 'dx' + str(dx)+ 'H1' +str(H1)  +'Salt.bin'
    writefield(file_name, salt)

def write_eta(dx,H1,eta):
    file_name = 'dx' + str(dx)+ 'H1' +str(H1)  +'Eta.bin'
    writefield(file_name, eta)

def write_uvel(dx,H1,uvel):
    file_name = 'dx' + str(dx) +'H1' +str(H1) +'U.bin'
    writefield(file_name, uvel)

def write_vvel(dx,H1,vvel):
    file_name = 'dx' + str(dx) +'H1' +str(H1) +'V.bin'
    writefield(file_name, vvel)

def write_taux(dx,taux):
    file_name = 'dx' + str(dx) +'taux.bin'
    writefield(file_name, taux)

def write_tauy(dx,tauy):
    file_name = 'dx' + str(dx) +'tauy.bin'
    writefield(file_name, tauy)
