import numpy as np
import datetime
import shutil
#Import the other module in this directory
import arcf_initialconds as arcf

#Define the vertical grid
#Set the the level thicknesses
z1 = np.concatenate([-0.5*np.ones((50)) ,np.array([-0.6, -0.7, -0.8, -0.9]),
                   -1*np.ones((14)),-1.5*np.ones((32))] )
#Take the cumulative sum to get the depths
z=np.cumsum(z1)

#Set the parameters
dx = 5e2 #Zonal grid spacing metres
dy = 5e2 #Meridional grid spacing metres
nx = 128 #Number of zonal grid points
ny = 768 #Number of meridional grid points
f = 1.4e-4 #Coriolis parameter
H1 = 30 #Stratification depth
dh = 5 #Stratification thickness
dSh = 1.25 #Salinity difference across the front
Lf = 1e4 #Frontal width metres

#Create the input files
arcf.create_arcf_binary(z,dx, dy, nx, ny, f, H1, dh, dSh, Lf)

#Create the datestamp
date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
time_stamp = date[0:4] + '_' + date[5:7] + '_'+ date[8:10] + '_'+ date[11:13] + '_'+ date[14:16]

#Create a datestamped version of the python scripts used to generate the inputs
## YOU need to adjust the path for this
shutil.copy2('/path/to/arcf_initialconds.py', 'arcf_initialconds_' + time_stamp + '.py')
