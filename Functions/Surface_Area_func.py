#!/home/jcibanezm/codes/ytt30/yt-x86_64/bin/python

from yt.config import ytcfg;ytcfg["yt","__withinreason"]="True"
import os
import uuid

# # Compute surface Area
# This notebook reads the cloud object information and recnstructs a 3D Object with the cloud information.
# Then it looks for the surface Pixels and computes the surface Area of the cloud.


## Convert the notebook to a python script this easy:
# ipython nbconvert --to python notebook.ipynb

# Calculate the accreted particles in the cloud
# This notebook takes a cloud, looks for the particles inside at each timestep and makes an array of particles in cloud as a function of time.
# It then looks backwards and asks which particles where accreted at each timestep and makes an array called accreted particles here, as a function of time.
# As soon as it reaches the beginning of the simulation, it takes all the particles in the cloud at the final step, checks all the accreted particles and maes a particle history for each particle. (gotta think about this step in more detail.)

import yt
import os
import math
from   yt import derived_field
import matplotlib.pyplot as plt
import numpy as np
from   matplotlib.colors import LogNorm
from   yt.units import pc, kpc, second, Kelvin, gram, erg, cm

import copy
from   matplotlib.colors import LogNorm

import timeit

# Define some constant parameters to be used.
mp      = 1.6726e-24  * gram # g
mu      = 1.2924
kb      = 1.3806e-16  *erg / Kelvin # erg K-1
GNewton = 6.6743e-8   * cm**3 / (gram * second**2 )# cm3 g-1 s-2
Msun    = 1.9884e33   * gram
mm      = mu*mp

ppc = 3.08567758e18

# Create a derived field.
@derived_field(name="numdens", units="1/cm**3", force_override=True)
def numdens(field, data):
    dens_here = data["dens"].value
    dens_here = dens_here * gram / cm**3
    return dens_here/mm

yt.add_field('numdens', function=numdens, units="1/cm**3", force_override=True)


def read_cloud_data(Cloud_name, resolution, snap, cloud_pop=0, SG=True):
    
    import numpy as np
    
    if cloud_pop == 1:
        cloud_data_dir = "/data/gamera/jcibanezm/StratBox/AccretionPaper/CloudPopulation/CloudObjects/" 
        if SG:
            # Navigate the folder names.
            snp_here = (int((int(Cloud_name.split("_")[1].split("t")[1]) + 2 ) / 10. ) - 230) * 10
            snp_str  = "%.2i" %snp_here
            sub_dir  = "t230_snpshot" + snp_str + "_ncut100"
            
            cloud_data_dir +="SG/09pc/" + sub_dir 
        else:
            # Navigate the folder names.
            snp_here = (int((int(Cloud_name.split("_")[1].split("t")[1]) + 2 ) / 10. ) - 230) * 10 + 10
            snp_str  = "%.2i" %snp_here
            print (snp_str)
            sub_dir  = "t230_snpshot" + snp_str + "_ncut100"

            #cloud_data_dir +="SG/09pc/" + sub_dir

            cloud_data_dir +="NoSG/1pc/" + sub_dir
            # Given the Cloud name, e.g. c000_t2300_n100, I should go to the proper folder and load the adequate file...
            # Test Do it and test it!!
        
        directory = cloud_data_dir
        filename  = Cloud_name
        
    else:
        # High resolution clouds case.
        cloud_data_dir = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Cloud_Object" 
        directory      = cloud_data_dir + "/%s" %Cloud_name + "/%.2ipc" %(resolution*10)
        filename       = Cloud_name + "_%.2ipc" %(resolution*10) + "_snp%.2i" %snap 

    # Directory where cloud properties are located.

    if not os.path.exists(directory):
        print "Not such directory, are you sure the cloud properties are stored there ?"
    
    # name of the cloud properties file.
    name = filename+".dat"

    # Basic variables.
    Header         = []
    num_cells      = 0
    num_properties = 0

    # Open the file and start reading line by line.
    f = open("%s/%s"%(directory,name), "r")
    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())
    
    # Blank space in this files only.
    if cloud_pop == 1: f.readline()

    # Read the number of clouds. (I could put a conditional that I had indeed red the number of clouds)
    entry1     = f.readline()

    cname = str(entry1.split("\t")[1].split("\n")[0])

    if cname != Cloud_name:
        print ("Something wrong here. I'm reading cloud %s, and was supposed to be reading cloud %s" %(cname, Cloud_name))

    # Read the number of properties.
    entry2     = f.readline().split("\t")
    num_properties = int(entry2[1].split("\n")[0])

    entry3 = f.readline().split("\t")
    num_cells = int(entry3[1].split("\n")[0])
    
    if cloud_pop == 1: f.readline()

    # Initialize the arrays of the names of the cloud properties and their values.
    properties_values      = []
    cloud_properties_names = ()

    # Read the names of the properties
    props_read = f.readline()

    props = (props_read.split()[0], props_read.split()[1], props_read.split()[2], props_read.split()[3], 
             props_read.split()[4], props_read.split()[5], props_read.split()[6], props_read.split()[7], 
             props_read.split()[8], props_read.split()[9], props_read.split()[10], props_read.split()[11], 
             props_read.split()[12], props_read.split()[13], props_read.split()[14], props_read.split()[15])

    # Make the properties arrays
    for i in range(num_properties):
        exec("%s = []" %(props[i]))

    # Read the values of each property and save it in the appropriate array.
    for i in range(num_cells):
        values = f.readline()
        for j in range(num_properties):
            exec("%s.append(float(%e))" %(props[j], float(values.split()[j])))

    for i in range(num_properties):
        exec("%s = np.array(%s)" %(props[i], props[i]))
        
    cloud_data = {'info':'Cloud data here.'}
    
    for this_prop in props:
        exec("cloud_data['%s'] = %s " %(this_prop, this_prop))
        
    if cloud_pop == 1: 
        if 'numdens' not in props:
            numdens = cloud_data["dens"] / mm.value
            cloud_data["numdens"] = numdens
    
    return cloud_data



def Calculate_Cloud_Surface_Area(Cloud_name, resolution, snapshot, cloud_pop=0, SG=True):
    """
    Given a cloud name, the resolution and the snapshot, This function restores the cloud data and computes the surface Area.
    Works for AMR grids.
    """
    
    Cloud_data = read_cloud_data(Cloud_name, resolution, snapshot, cloud_pop=cloud_pop, SG=SG)
    
    # Reconstruct the Grid in a Uniform Grid Structure
    dx   = np.min(Cloud_data["dx"])
    
    xmin = np.min(Cloud_data["x"]) - 2*dx
    xmax = np.max(Cloud_data["x"]) + 2*dx
    ymin = np.min(Cloud_data["y"]) - 2*dx
    ymax = np.max(Cloud_data["y"]) + 2*dx
    zmin = np.min(Cloud_data["z"]) - 2*dx
    zmax = np.max(Cloud_data["z"]) + 2*dx
    
    Lx = xmax - xmin 
    Ly = ymax - ymin
    Lz = zmax - zmin
    
    Nx = int(Lx / dx + 2)
    Ny = int(Ly / dx + 2)
    Nz = int(Lz / dx + 2)
    
    UG = np.zeros((Nx, Ny, Nz))
    
    print ("===============================================")
    print ("Reconstructing a uniform grid for this cloud")
    print ("Uniform Grid properties:")
    print ("Lz = %.2f   Ly = %.2f   Lz = %.2f"%(Lx/ppc, Ly/ppc, Lz/ppc))
    print ("dx = %.2f pc" %(dx/ppc))
    print ("Nx = %i,    Ny = %i,    Nz = %i"  %(Nx, Ny, Nz))
    print ("xmin = %.2f      xmax = %.2f"     %(xmin/ppc, xmax/ppc))
    print ("ymin = %.2f      ymax = %.2f"     %(ymin/ppc, ymax/ppc))
    print ("zmin = %.2f      zmax = %.2f"     %(zmin/ppc, zmax/ppc))
    print ("===============================================")
    
    
    # Initialize a dictionary that will contain the grid where the cloud will be mapped.
    Cloud_Grid = {"hello": "I'm the new grid with the cloud information"}
    	
    # Initialize basic properties of the uniform grid for the cloud. Reconstruct the x, y and z positions of the cloud uniform grid.
    # Also inizialize the cell size dx, dy, dz.	
    X  = np.zeros_like(UG)
    DX = np.zeros_like(UG)	
    for i in range(np.shape(UG)[0]):
        X[i,:,:]  = xmin + i*dx
        DX[i,:,:] = dx
    
    Cloud_Grid["x"]  = X
    Cloud_Grid["dx"] = DX	
    del X
    
    Y = np.zeros_like(UG)	
    for j in range(np.shape(UG)[1]):
        Y[:,j,:]  = ymin + i*dx
        DX[:,j,:] = dx
        
    Cloud_Grid["y"]  = Y
    Cloud_Grid["dy"] = DX	
    del Y
    
    Z = np.zeros_like(UG)	
    for k in range(np.shape(UG)[2]):
        Z[:,:,k]  = zmin + i*dx
        DX[:,:,k] = dx
        
    Cloud_Grid["z"] = Z
    Cloud_Grid["dz"] = DX	
    del Z
    del DX
    
    
    # Read a position of the cloud, get the corresponding index in the uniform grid and save the appropriate information.
    ndens = np.zeros_like(UG)
    magx  = np.zeros_like(UG)
    magy  = np.zeros_like(UG)
    magz  = np.zeros_like(UG)
    velx  = np.zeros_like(UG)
    vely  = np.zeros_like(UG)
    velz  = np.zeros_like(UG)
    temp  = np.zeros_like(UG)
    
    for cc in range(len(Cloud_data["x"])):
        
        DX2dx = int(Cloud_data["dx"][cc]/dx)
        
        if DX2dx == 1:
            iindex   = int((Cloud_data["x"][cc] - xmin ) / dx + 0.1)
            jindex   = int((Cloud_data["y"][cc] - ymin ) / dx + 0.1 )
            kindex   = int((Cloud_data["z"][cc] - zmin ) / dx + 0.1 )
            
            ndens[iindex, jindex, kindex] = Cloud_data["numdens"][cc]
            
        else:
            # I have to cycle Over the AMR grid.
            for sub_cycle in range(DX2dx**3):
                # Convert from a 1D array to a 3D array indexes.
                kk =  sub_cycle / (DX2dx * DX2dx)
                jj = (sub_cycle - kk*DX2dx*DX2dx) / DX2dx
                ii =  sub_cycle - kk*DX2dx*DX2dx - jj*DX2dx
                
                # Get the proper index with respect to the box for the sub-cycling amr grid.
                #iindex   = int((Cloud_data["x"][cc] - xmin - Cloud_data["dx"][cc]/2.0 + (ii*2+1)*dx/2.) / dx + 0.1 )
                #jindex   = int((Cloud_data["y"][cc] - ymin - Cloud_data["dy"][cc]/2.0 + (jj*2+1)*dx/2.) / dx + 0.1 )
                #kindex   = int((Cloud_data["z"][cc] - zmin - Cloud_data["dz"][cc]/2.0 + (kk*2+1)*dx/2.) / dx + 0.1 )
    
                iindex   = int((Cloud_data["x"][cc] - xmin - Cloud_data["dx"][cc]/2.0 + dx/2.) / dx + 0.1 ) + ii
                jindex   = int((Cloud_data["y"][cc] - ymin - Cloud_data["dy"][cc]/2.0 + dx/2.) / dx + 0.1 ) + jj
                kindex   = int((Cloud_data["z"][cc] - zmin - Cloud_data["dz"][cc]/2.0 + dx/2.) / dx + 0.1 ) + kk
                
                # Store the local density value here.
                ndens[iindex, jindex, kindex] = Cloud_data["numdens"][cc]
                temp[iindex, jindex, kindex]  = Cloud_data["temp"][cc]
                
                velx[iindex, jindex, kindex]  = Cloud_data["velx"][cc]
                vely[iindex, jindex, kindex]  = Cloud_data["vely"][cc]
                velz[iindex, jindex, kindex]  = Cloud_data["velz"][cc]
    
                magx[iindex, jindex, kindex]  = Cloud_data["magx"][cc]
                magy[iindex, jindex, kindex]  = Cloud_data["magy"][cc]
                magz[iindex, jindex, kindex]  = Cloud_data["magz"][cc]
    
    
    Cloud_Grid["numdens"] = ndens
    Cloud_Grid["temp"]    = temp
    
    Cloud_Grid["magx"] = magx
    Cloud_Grid["magy"] = magy
    Cloud_Grid["magz"] = magz
    Cloud_Grid["velx"] = velx
    Cloud_Grid["vely"] = vely
    Cloud_Grid["velz"] = velz
    
    del ndens 
    del magx  
    del magy  
    del magz  
    del velx  
    del vely  
    del velz
    del temp
    
    # Make an array asking if I am a surface pixel or not.
    # Compute how many faces of this pixel correspond to a surface face.
    Surface_pixel = np.zeros_like(UG)
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                if Cloud_Grid["numdens"][i, j, k] != 0:
                    # I'm in a cloud pixel.
                    # Now check my neighbours and see if I'm inside the cloud or a surface pixel.
                    if Cloud_Grid["numdens"][i+1, j, k] == 0 : Surface_pixel[i,j,k] +=1
                    if Cloud_Grid["numdens"][i-1, j, k] == 0 : Surface_pixel[i,j,k] +=1
                    if Cloud_Grid["numdens"][i, j+1, k] == 0 : Surface_pixel[i,j,k] +=1
                    if Cloud_Grid["numdens"][i, j-1, k] == 0 : Surface_pixel[i,j,k] +=1
                    if Cloud_Grid["numdens"][i, j, k+1] == 0 : Surface_pixel[i,j,k] +=1
                    if Cloud_Grid["numdens"][i, j, k-1] == 0 : Surface_pixel[i,j,k] +=1
    
    Cloud_Grid["surface_pixels"] = Surface_pixel
    
    # Compute the total surface Area
    total_surface_area = np.sum(Cloud_Grid["surface_pixels"]*Cloud_Grid["dx"]**2)
    
    # Store this value in the Cloud object.
    Cloud_Grid["total_surface"] = total_surface_area
    
    return Cloud_Grid
