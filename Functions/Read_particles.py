#!/home/jcibanezm/codes/ytt30/yt-x86_64/bin/python

from yt.config import ytcfg;ytcfg["yt","__withinreason"]="True"
import os
import uuid

## Convert the notebook to a python script this easy:
# ipython nbconvert --to python notebook.ipynb


# This notebook is supposed to be converted into a function that will read the particles inside each cloud at each timestep, given the name of the cloud, the resolution nad the start and end time of the simulation and return a bunch of arrays.

import yt
import os
import math
from yt import derived_field
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from yt.units import pc, kpc, second, Kelvin, gram, erg, cm

import copy
from matplotlib.colors import LogNorm

import timeit


# Define some constant parameters to be used.
mp      = 1.6726e-24  * gram # g
mu      = 1.2924
kb      = 1.3806e-16  *erg / Kelvin # erg K-1
GNewton = 6.6743e-8   * cm**3 / (gram * second**2 )# cm3 g-1 s-2
Msun    = 1.9884e33   * gram
mm      = mu*mp

ppc     = 3.08e18

# Create a derived field.
@derived_field(name="numdens", units="1/cm**3", force_override=True)
def numdens(field, data):
    dens_here = data["dens"].value
    dens_here = dens_here * gram / cm**3
    return dens_here/mm

yt.add_field('numdens', function=numdens, units="1/cm**3", force_override=True)

####
# Read the information here
####
def Read_Particles_snapshot(Cloud_name, resolution, tnow, ncut=100):
    """
    Read the particles inside the cloud.
    """
    
    print("Reading the particle information of cloud %s, with resolution %1.2f pc, at time %.2f Myr and density threshold %i cc" %(Cloud_name, resolution, tnow, ncut))

    if resolution < 0.1:
        resolution_str = "%.3i" %(resolution*100)
    else:
        resolution_str = "%.2i" %(resolution*10)
    
    save_particles_dir = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Particles/%s/%spc/%incut" %(Cloud_name, resolution_str, ncut) 

    if not os.path.exists(save_particles_dir):
        print ("The directory you are looking for doesn't exist.")

    snapshot     = int(tnow*10)
    snapshot_str = str(snapshot)

    if (len(snapshot_str) == 1):
        snapshot_str = "0"+ snapshot_str
    fn2 = "snpsht"+snapshot_str

    # concatenate the filename.
    filename = Cloud_name + "_PartInCloud_" + fn2

    name = filename+".dat"
    f = open("%s/%s"%(save_particles_dir,name), 'r')

    # Read the header of the file
    Header = []

    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())

    # Read the number of particles
    entry0 = f.readline()
    num_particles = int(entry0.split()[1])
    
    # Read the number of properties of the particles
    entry1 = f.readline()
    num_properties = int(entry1.split()[1])
    
    # Read the names of the properties
    props_read = f.readline()
    
    props = []    
    for i in range(num_properties):
        props.append(props_read.split()[i*2])
        
    # Make the properties arrays
    #for i in range(num_properties):
    #    exec("%s = []" %(props[i]))
    
    # Make the properties arrays
    for i in range(num_properties):
        exec("%s = []" %(props[i]))
    
    # Read the values of each property and save it in the appropriate array.
    for i in range(num_particles):
        values = f.readline()
        for j in range(num_properties):
            if props[j] == "tag":
                exec("%s.append(int(%i))" %(props[j], int(values.split()[j])))
            else:
                exec("%s.append(float(%e))" %(props[j], float(values.split()[j])))
    
    f.close()
    
    # Make a dictionary to return the data.
    for j in range(num_properties):
        if j == 0 :
            exec("particle_data = {'%s': %s} " %(props[j], props[j]))
        else:
            exec("particle_data['%s'] = %s" %(props[j], props[j]))
    
    return particle_data


# ---
# # Make an array of arrays
def get_particles_evolution(Cloud_name, resolution, tinit, tfinal, ncut=100):
    
    global_particles = {'info':'This dictionary contains the information of the particles inside the cloud for all snapshots.  in order to get the key names type global_particles.keys(). the first entry corresponds to the snapshot (time) and the  second corresponds to the particle entries.'}

    global_particles["time"] = []
    
    particle_data = Read_Particles_snapshot(Cloud_name, resolution, 0., ncut)
    for key_name in particle_data.keys():
        global_particles[key_name] = []

    # Number of snapshots.
    nsnaps = int((tfinal - tinit)*10)
    
    for i in range(nsnaps+1):
        tnow = tinit + i/10.
        particle_data = Read_Particles_snapshot(Cloud_name, resolution, tnow, ncut)
        global_particles["time"].append(tnow)
        for key_name in particle_data.keys():
            global_particles[key_name].append(particle_data['%s'%key_name])

    return global_particles


############################
# Read the information here
############################

def Read_Inst_Acc_Particles_snapshot(Cloud_name, resolution, tnow, ncut=100):
    """
    Read the information of the particles at two snapshots, compute the accreted particles between snapshots and calculate the properties of the accreted particles.
    Input: Cloud_name, resolution, evolutionary time, and density threshold.
    """
    
    print("Reading the properties of the particles accreted between two subsequent snapshots. cloud %s, with resolution %1.2f pc, at time %.2f Myr and density threshold %i cc" %(Cloud_name, resolution, tnow, ncut))


    if resolution < 0.1:
        resolution_str = "%.3i" %(resolution*100)
    else:
        resolution_str = "%.2i" %(resolution*10)

    save_particles_dir = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Inst_Acc_Particles/%s/%spc/%incut" %(Cloud_name, resolution_str, ncut) 

    if not os.path.exists(save_particles_dir):
        print ("The directory you are looking for doesn't exist.")

    snapshot     = int(tnow*10)
    snapshot_str = str(snapshot)

    if (len(snapshot_str) == 1):
        snapshot_str = "0"+ snapshot_str
    fn2 = "snpsht"+snapshot_str

    # concatenate the filename.
    filename = Cloud_name + "_InstAccParts_" + fn2
    
    name = filename+".dat"
    f = open("%s/%s"%(save_particles_dir,name), 'r')

    # Read the header of the file
    Header = []

    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())

    # Read the number of particles
    entry0 = f.readline()
    num_particles = int(entry0.split()[1])

    # Read the number of properties of the particles
    entry1 = f.readline()
    num_properties = int(entry1.split()[1])

    # Read the names of the properties
    props_read = f.readline()

    props = []
    for i in range(num_properties):
        props.append(props_read.split()[i*2])

    # Make the properties arrays
    #for i in range(num_properties):
    #    exec("%s = []" %(props[i]))

    # Make the properties arrays
    for i in range(num_properties):
        exec("%s = []" %(props[i]))

    # Read the values of each property and save it in the appropriate array.
    for i in range(num_particles):
        values = f.readline()
        for j in range(num_properties):
            if props[j] == "tag":
                exec("%s.append(int(%i))" %(props[j], int(values.split()[j])))
            else:
                exec("%s.append(float(%e))" %(props[j], float(values.split()[j])))

    f.close()

    # Make a dictionary to return the data.
    for j in range(num_properties):
        if j == 0 :
            exec("particle_data = {'%s': %s} " %(props[j], props[j]))
        else:
            exec("particle_data['%s'] = %s" %(props[j], props[j]))


    return particle_data


def get_Inst_Acc_particles_evolution(Cloud_name, resolution, tinit, tfinal, ncut=100):
    
    global_particles = {'info':'This dictionary contains the information of the particles inside the cloud for all snapshots. \
 in order to get the key names type global_particles.keys(). the first entry corresponds to the snapshot (time) and the \
 second corresponds to the particle entries.'}

    global_particles["time"] = []
    
    particle_data = Read_Inst_Acc_Particles_snapshot(Cloud_name, resolution, 0., ncut)
    for key_name in particle_data.keys():
        global_particles[key_name] = []

    # Number of snapshots.
    nsnaps = int((tfinal - tinit)*10)
    
    for i in range(nsnaps+1):
        tnow = tinit + i/10.
        particle_data = Read_Inst_Acc_Particles_snapshot(Cloud_name, resolution, tnow, ncut)
        global_particles["time"].append(tnow)
        for key_name in particle_data.keys():
            global_particles[key_name].append(particle_data['%s'%key_name])

    return global_particles


#####################################################################################
# Read the global accretion particles.
#####################################################################################

############################
# Read the information here
############################

def Read_Acc_Particles_snapshot(Cloud_name, resolution, tnow, ncut=100):
    """
    Reading the total particles accreted onto the cloud.
    """
    
    print("Reading particle properties of all accreted particles. cloud %s, with resolution %1.2f pc, at time %.2f Myr and density threshold %i cc" %(Cloud_name, resolution, tnow, ncut))


 
    if resolution < 0.1:
        resolution_str = "%.3i" %(resolution*100)
    else:
        resolution_str = "%.2i" %(resolution*10)

    save_particles_dir = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Acc_Particles/%s/%spc/%incut" %(Cloud_name, resolution_str, ncut) 

    if not os.path.exists(save_particles_dir):
        print ("The directory %s doesn't exist." %save_particles_dir)

    snapshot     = int(tnow*10)
    snapshot_str = str(snapshot)

    if (len(snapshot_str) == 1):
        snapshot_str = "0"+ snapshot_str
    fn2 = "snpsht"+snapshot_str

    # concatenate the filename.
    filename = Cloud_name + "_AccParts_" + fn2
    
    name = filename+".dat"
    f = open("%s/%s"%(save_particles_dir,name), 'r')

    # Read the header of the file
    Header = []

    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())

    # Read the number of particles
    entry0 = f.readline()
    num_particles = int(entry0.split()[1])

    # Read the number of properties of the particles
    entry1 = f.readline()
    num_properties = int(entry1.split()[1])

    # Read the names of the properties
    props_read = f.readline()

    props = (props_read.split()[0], props_read.split()[1], props_read.split()[3], props_read.split()[5], 
             props_read.split()[7], props_read.split()[9], props_read.split()[11], props_read.split()[13], 
             props_read.split()[15], props_read.split()[17], props_read.split()[19], props_read.split()[21])
    
    # Make the properties arrays
    for i in range(num_properties):
        exec("%s = []" %(props[i]))

    # Read the values of each property and save it in the appropriate array.
    for i in range(num_particles):
        values = f.readline()
        for j in range(num_properties):
            if props[j] == "PinC_tag":
                exec("%s.append(int(%i))" %(props[j], int(values.split()[j])))
            else:
                exec("%s.append(float(%f))" %(props[j], float(values.split()[j])))

    f.close()

    for j in range(num_properties):
        if j == 0 :
            exec("particle_data = {'%s': %s} " %(props[j], props[j]))
        else:
            exec("particle_data['%s'] = %s" %(props[j], props[j]))
    
    return particle_data


def get_Acc_particles_evolution(Cloud_name, resolution, tinit, tfinal, ncut=100):
    
    global_particles = {'info':'This dictionary contains the information of the particles inside the cloud for all snapshots. \
 in order to get the key names type global_particles.keys(). the first entry corresponds to the snapshot (time) and the \
 second corresponds to the particle entries.'}

    global_particles["time"] = []
    
    particle_data = Read_Acc_Particles_snapshot(Cloud_name, resolution, 0.1, ncut)
    for key_name in particle_data.keys():
        global_particles[key_name] = []

    # Number of snapshots.
    nsnaps = int((tfinal - tinit)*10)
    
    for i in range(nsnaps+1):
        
        tnow          = tinit + i/10.
        particle_data = Read_Acc_Particles_snapshot(Cloud_name, resolution, tnow, ncut)
        global_particles["time"].append(tnow)
        for key_name in particle_data.keys():
            global_particles[key_name].append(particle_data['%s'%key_name])

    return global_particles


##################################################################
# get the cloud properties here.

def get_cloud_data(Cloud_name, resolution, snap, nmin, nmax):
    print("Reading the cloud properties. Cloud %s, resolution %.2f pc, at time %.2f, between thresholds %.2f and %.2f cm-3" %(Cloud_name, resolution, snap/10., nmin, nmax))
    cloud_data_dir = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Cloud_Object"

    directory = cloud_data_dir + "/%s" %Cloud_name + "/%.2ipc" %(resolution*10)

    # Directory where cloud properties are located.

    if not os.path.exists(directory):
        print "Not such directory, are you sure the cloud properties are stored there ?"

    filename = Cloud_name + "_%.2ipc" %(resolution*10) + "_snp%.2i" %snap 

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
            exec("%s.append(float(%f))" %(props[j], float(values.split()[j])))

    for i in range(num_properties):
        exec("%s = np.array(%s)" %(props[i], props[i]))

    ##################################################
    # Restore the cloud properties for a density range.
    ##################################################

    within_density_range = []
    for j in range(len(numdens)):
        if ( numdens[j] >= nmin and numdens[j] <= nmax):
            within_density_range.append(1)
        else:
            within_density_range.append(0)

    within_density_range = np.array(within_density_range)

    mass = sum(cell_mass * within_density_range)

    if mass ==  0:
        mass = 1.0e-99

    volume        = sum(cell_volume * within_density_range)
    cell_num      = sum(within_density_range) 
    avg_dens      = mass / volume
    spherical_rad = (3 * volume / (4 * math.pi))**(1.0/3.0)
    tff           = math.sqrt( 3*math.pi / ( 32 * GNewton * avg_dens ) )
    CM_x          = sum(x * cell_mass * within_density_range) / mass
    CM_y          = sum(y * cell_mass * within_density_range) / mass
    CM_z          = sum(z * cell_mass * within_density_range) / mass

    x_min         = min(x * within_density_range)
    y_min         = min(y * within_density_range)
    z_min         = min(z * within_density_range)

    x_max         = max(x * within_density_range)
    y_max         = max(y * within_density_range)
    z_max         = max(z * within_density_range)

    bulk_vel_x    = sum(velx* cell_mass * within_density_range) / mass
    bulk_vel_y    = sum(vely* cell_mass * within_density_range) / mass
    bulk_vel_z    = sum(velz* cell_mass * within_density_range) / mass

    velocity_magnitude = np.sqrt(velx**2 + vely**2 + velz**2)
    velocity_magnitude = np.array(velocity_magnitude)

    bulk_vel_3D   = sum(velocity_magnitude * cell_mass * within_density_range) / mass

    x_vel_disp    = np.sqrt(sum((velx - bulk_vel_x)**2 * cell_mass * within_density_range) / mass )
    y_vel_disp    = np.sqrt(sum((vely - bulk_vel_y)**2 * cell_mass * within_density_range) / mass )
    z_vel_disp    = np.sqrt(sum((velz - bulk_vel_z)**2 * cell_mass * within_density_range) / mass )

    vel_disp_3D   = math.sqrt(sum((velocity_magnitude - bulk_vel_3D)**2 * 
                                   cell_mass * within_density_range) / mass )

    avg_cs = sum(np.sqrt(temp * kb.value * 5.0 / 3.0 / mm.value ) * cell_mass * within_density_range ) / mass

    vel_disp_total = math.sqrt(avg_cs**2 + 1.0/3.0 * vel_disp_3D**2)

    alpha_virial = 5 * vel_disp_total**2 * spherical_rad / ( GNewton.value * mass)


    new_props = ('mass', 'volume', 'cell_num', 'avg_dens', 'spherical_rad', 'tff', 'CM_x', 'CM_y', 'CM_z',
                 'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max', 'bulk_vel_x', 'bulk_vel_y', 'bulk_vel_z',
                 'bulk_vel_3D', 'x_vel_disp', 'y_vel_disp', 'z_vel_disp', 'vel_disp_3D', 'avg_cs', 'vel_disp_total',
                 'alpha_virial')

    cloud_props = {'info':'Cloud properties here.'}

    for this_prop in new_props:
        exec("cloud_props['%s'] = %s " %(this_prop, this_prop))

    return cloud_props


def get_cloud_global_evolution(Cloud_name, resolution, tinit, tfinal, nmin, nmax):
    """
    Get the cloud properties between tinit and tfinal.
    """ 
    cloud_global_evolution = {'info': 'Global evolution of the cloud properties.'}
    
    cloud_global_evolution["time"] = []
    
    cloud_here = get_cloud_data(Cloud_name, resolution, 0., nmin, nmax)
    for key_name in cloud_here.keys():
        cloud_global_evolution[key_name] = []
    
    print ("Done reading and recovering the cloud properties. The available keys are:")
    print cloud_global_evolution.keys()

    nsnaps = int((tfinal - tinit)*10) 
    for snap in range(nsnaps+1):
        cloud_here = get_cloud_data(Cloud_name, resolution, snap, nmin, nmax)
        cloud_global_evolution["time"].append(snap/10.0)
        for key_name in cloud_here.keys():
            cloud_global_evolution[key_name].append(cloud_here['%s'%key_name])
        
    return cloud_global_evolution
