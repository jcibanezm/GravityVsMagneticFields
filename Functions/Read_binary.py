############################################################################
# Read the cloud data
############################################################################

# Read one snapshot

def Read_One_Snapshot(Cloud_name, resolution, snapshot, dir='default'):
    """
    This function reads the binary where the cloud object information is stored.
    """
    
    import cPickle as cP
    
    # Set the string for the resolution.
    if resolution < 0.1:
        resolution_str = "%.3i" %(resolution*100)
    else:
        resolution_str = "%.2i" %(resolution*10)
    
    if dir == 'default':
        # Set the directory where the cloud data is saved
        saved_cloud = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Cloud_Object/Binary/"
        saved_cloud = "%s/%s/%spc" %(saved_cloud, Cloud_name, resolution_str)
    else:
        saved_cloud = dir    

    print("Reading the File %s_CloudObject_snp%.3i.dat" %(Cloud_name, snapshot) )
    
    # Open the binary file and read the information within.
    p = open("%s/%s_CloudObject_snp%.3i.dat" %(saved_cloud, Cloud_name, snapshot), 'rb')
    Cloud_data = cP.load(p)
    p.close()
    
    # Return the information.
    return Cloud_data


def Read_Box_Snapshot(Cloud_name, resolution, snapshot, dir='default'):
    """
    This function reads the binary where the cloud object information is stored.
    """

    import cPickle as cP

    # Set the string for the resolution.
    if resolution < 0.1:
        resolution_str = "%.3i" %(resolution*100)
    else:
        resolution_str = "%.2i" %(resolution*10)

    if dir == 'default':
        # Set the directory where the cloud data is saved
        saved_cloud = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Cloud_Object/Binary/"
        saved_cloud = "%s/%s/%spc" %(saved_cloud, Cloud_name, resolution_str)
    else:
        saved_cloud = dir

    print("Reading the File %s_Box_snp%.3i.dat" %(Cloud_name, snapshot) )

    # Open the binary file and read the information within.
    p = open("%s/%s_Box_snp%.3i.dat" %(saved_cloud, Cloud_name, snapshot), 'rb')
    Cloud_data = cP.load(p)
    p.close()

    # Return the information.
    return Cloud_data


#####################################################
# Make an array of arrays for the cloud evolution
#####################################################

def Read_Cloud_Evolution(Cloud_name, resolution, tinit, tfin):
    
    global_data = {'info':'This dictionary contains the information of the cloud data as and its time evolution.'}

    global_data["time"] = []

    snapshot_data0 = Read_One_Snapshot(Cloud_name, resolution, 0.)
    for key_name in snapshot_data0.keys():
        global_data[key_name] = []

    # Number of snapshots.
    nsnaps = int((tfin - tinit)*10)

    for i in range(nsnaps+1):
        tnow = tinit + i/10.
        snapshot_data = Read_One_Snapshot(Cloud_name, resolution, i)
        global_data["time"].append(tnow)
        for key_name in snapshot_data0.keys():
            global_data[key_name].append(snapshot_data['%s'%key_name])
            
    return global_data





#################################################################################################################
# Read the particles in cloud
##################################################################################################################

def Read_One_Particle_Snapshot(Cloud_name, resolution, snapshot):
    """
    This function reads the binary where the cloud object information is stored.
    """
    
    import cPickle as cP
    
    # Set the string for the resolution.
    if resolution < 0.1:
        resolution_str = "%.3i" %(resolution*100)
    else:
        resolution_str = "%.2i" %(resolution*10)
    
    # Set the directory where the cloud data is saved
    
    saved_particles = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Particles/Binary"
    saved_particles = "%s/%s/%spc" %(saved_particles, Cloud_name, resolution_str)
    
    print("Reading the File %s_Particles_snp%.3i.dat" %(Cloud_name, snapshot) )
    
    # Open the binary file and read the information within.
    p = open("%s/%s_Particles_snp%.3i.dat" %(saved_particles, Cloud_name, snapshot), 'rb')
    particles_data = cP.load(p)
    p.close()
    
    # Return the information.
    return particles_data


########################################################
# Make an array of arrays for the particles evolution
########################################################

def Read_Particles_Evolution(Cloud_name, resolution, tinit, tfin):
    
    global_data = {'info':'This dictionary contains the information of the cloud data as and its time evolution.'}

    global_data["time"] = []

    snapshot_data0 = Read_One_Particle_Snapshot(Cloud_name, resolution, 0)
    for key_name in snapshot_data0.keys():
        global_data[key_name] = []

    # Number of snapshots.
    nsnaps = int((tfin - tinit)*10)

    for i in range(nsnaps+1):
        tnow = tinit + i/10.
        snapshot_data = Read_One_Particle_Snapshot(Cloud_name, resolution, i)
        global_data["time"].append(tnow)
        for key_name in snapshot_data0.keys():
            global_data[key_name].append(snapshot_data['%s'%key_name])
            
    return global_data


################################################################
# Read one snapshot and restore the cloud properties
################################################################
def Restore_One_Snapshot(Cloud_name, resolution, snapshot, nmin=100, nmax=1.0e99, object_here='clump'):
    """
    This function reads the binary where the cloud object information is stored.
    """
    
    import cPickle as cP
    import math
    import numpy as np
    
    # Set the string for the resolution.
    if resolution < 0.1:
        resolution_str = "%.3i" %(resolution*100)
    else:
        resolution_str = "%.2i" %(resolution*10)
    
    # Set the directory where the cloud data is saved
    saved_cloud = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Cloud_Object/Binary/"
    saved_cloud = "%s/%s/%spc" %(saved_cloud, Cloud_name, resolution_str)
    
    print("Reading the File %s_CloudObject_snp%.3i.dat" %(Cloud_name, snapshot) )
    
    # Open the binary file and read the information within.
    p = open("%s/%s_CloudObject_snp%.3i.dat" %(saved_cloud, Cloud_name, snapshot), 'rb')
    Cloud_data = cP.load(p)
    p.close()
    
    
    ### Regenerate the cloud properties from the data in the cloud.
    # Define some constant parameters to be used.    
    GNewton = 6.67259e-8
    kb      = 1.380658e-16   
    mp      = 1.6726e-24  # g
    mu      = 1.2924
    mm      = mu*mp

    within_density_range = []
    above_density_threshold = []
    
    num_cells = len(Cloud_data["numdens"])
    for j in range(num_cells):
        if (Cloud_data["numdens"][j] >= nmin) :
            above_density_threshold.append(1)
            if Cloud_data["numdens"][j] <= nmax :
                within_density_range.append(1)
            else:
                within_density_range.append(0)
        else:
            above_density_threshold.append(0)
            within_density_range.append(0)    
    
    within_density_range    = np.array(within_density_range)
    above_density_threshold = np.array(above_density_threshold)        

    mass          = np.sum(Cloud_data["cell_mass"] * within_density_range)
    mass_total    = np.sum(Cloud_data["cell_mass"] )

    if mass ==  0:
        mass = 1.0e-99
        
    if object_here == 'clump':
        volume    = np.sum(Cloud_data["cell_volume"] * above_density_threshold)
    elif object_here == 'cloud':
        volume    = np.sum(Cloud_data["cell_volume"])

    volume_total  = np.sum(Cloud_data["cell_volume"])        
 
    cell_num      = np.sum(within_density_range) #len(clouds_restored[ii].dens)
    avg_dens      = mass_total / volume_total
    spherical_rad = (3. * volume / (4. * math.pi))**(1.0/3.0)
    tff           = np.sqrt( 3. * math.pi / ( 32. * GNewton * avg_dens ) )
    CM_x          = np.sum(Cloud_data["x"] * Cloud_data["cell_mass"] ) / mass
    CM_y          = np.sum(Cloud_data["y"] * Cloud_data["cell_mass"] ) / mass
    CM_z          = np.sum(Cloud_data["z"] * Cloud_data["cell_mass"] ) / mass
        
    x_min         = np.min(Cloud_data["x"] * within_density_range)
    y_min         = np.min(Cloud_data["y"] * within_density_range)
    z_min         = np.min(Cloud_data["z"] * within_density_range)
        
    x_max         = np.max(Cloud_data["x"] * within_density_range)
    y_max         = np.max(Cloud_data["y"] * within_density_range)
    z_max         = np.max(Cloud_data["z"] * within_density_range)
        
    bulk_vel_x    = np.sum(Cloud_data["velx"] * Cloud_data["cell_mass"] * within_density_range) / mass
    bulk_vel_y    = np.sum(Cloud_data["vely"] * Cloud_data["cell_mass"] * within_density_range) / mass
    bulk_vel_z    = np.sum(Cloud_data["velz"] * Cloud_data["cell_mass"] * within_density_range) / mass
    
    velocity_magnitude = np.zeros_like(Cloud_data["velx"])
    velocity_magnitude = np.sqrt(Cloud_data["velx"]*Cloud_data["velx"] + 
                                 Cloud_data["vely"]*Cloud_data["vely"] + 
                                 Cloud_data["velz"]*Cloud_data["velz"] )
    
    bulk_vel_3D   = np.sum(velocity_magnitude * Cloud_data["cell_mass"] * within_density_range) / mass
        
    x_vel_disp    = np.sqrt(np.sum((Cloud_data["velx"] - bulk_vel_x)**2 * Cloud_data["cell_mass"] * within_density_range) / mass )
    y_vel_disp    = np.sqrt(np.sum((Cloud_data["vely"] - bulk_vel_y)**2 * Cloud_data["cell_mass"] * within_density_range) / mass )
    z_vel_disp    = np.sqrt(np.sum((Cloud_data["velz"] - bulk_vel_z)**2 * Cloud_data["cell_mass"] * within_density_range) / mass )
        
    vel_disp_3D   = np.sqrt(np.sum((velocity_magnitude - bulk_vel_3D)**2 * Cloud_data["cell_mass"] * within_density_range) / mass )
       
    avg_cs = np.sum( np.sqrt(Cloud_data["temp"] *kb * 5.0 / 3.0 / mm ) * Cloud_data["cell_mass"] * within_density_range ) / mass
        
    vel_disp_total = np.sqrt(avg_cs**2 + 1.0/3.0 * vel_disp_3D**2)
         
    alpha_virial = 5 * vel_disp_total**2 * spherical_rad / ( GNewton * mass)
    
    Cloud_props = {"info": "This dictionary contains the properties of the restored cloud"}
    
    Cloud_props["mass"] = mass
    Cloud_props["volume"] = volume
    Cloud_props["CM_x"] = CM_x
    Cloud_props["CM_y"] = CM_y
    Cloud_props["CM_z"] = CM_z
    Cloud_props["BV_x"] = bulk_vel_x
    Cloud_props["BV_y"] = bulk_vel_y
    Cloud_props["BV_z"] = bulk_vel_z
    Cloud_props["BV_3D"] = bulk_vel_3D
    
    Cloud_props["vel_disp_3D"]    = vel_disp_3D
    Cloud_props["vel_disp_total"] = vel_disp_total
    Cloud_props["sound_speed"]    = avg_cs
    Cloud_props["alpha_virial"]   = alpha_virial

    Cloud_props["radius"] = spherical_rad
    Cloud_props["tff"]   = tff

    
    # compute the cloud properties as well !! Return 2 dictionaries.
    
    # Return the information.
    return Cloud_data, Cloud_props

def Restore_Cloud_Evolution(Cloud_name, resolution, tinit, tfin, nmin=100, nmax=1.0e99):
    
    global_data = {'info':'This dictionary contains the information of the cloud data as and its time evolution.'}
    global_props = {"info":"this dictionary contains the global properties as the cloud evolves"}
        
    global_data["time"]  = []
    global_props["time"] = []

    snapshot_data0, snapshot_props0 = Restore_One_Snapshot(Cloud_name, resolution, 0.)

    for key_name in snapshot_data0.keys():
        global_data[key_name] = []
        
    for key_name in snapshot_props0.keys():
        global_props[key_name] = []
        
    # Number of snapshots.
    nsnaps = int((tfin - tinit)*10)

    for i in range(nsnaps+1):
        tnow = tinit + i/10.
        
        snapshot_data, snapshot_props = Restore_One_Snapshot(Cloud_name, resolution, i, nmin=nmin, nmax=nmax)
        
        global_data["time"].append(tnow)
        global_props["time"].append(tnow)
        
        for key_name in snapshot_data0.keys():
            global_data[key_name].append(snapshot_data['%s'%key_name])
            
        for key_name in snapshot_props0.keys():
            global_props[key_name].append(snapshot_props['%s'%key_name])
    
    return global_data, global_props


#################################################################################################
#  Read instantaneous accretion.
#################################################################################################


def Read_One_InstAcc(Cloud_name, resolution, snapshot):
    """
    This function reads the binary where the cloud object information is stored.
    """
    
    import cPickle as cP
    
    # Set the string for the resolution.
    if resolution < 0.1:
        resolution_str = "%.3i" %(resolution*100)
    else:
        resolution_str = "%.2i" %(resolution*10)
    
    saved_inst_acc = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/Inst_Acc_Particles/%s/%spc/" %(Cloud_name, resolution_str,) 

    snapshot_str      = str(snapshot)

    if (len(snapshot_str) == 1):
        snapshot_str = "0"+ snapshot_str
    fn2 = "snpsht"+snapshot_str

    # concatenate the filename.
    filename = Cloud_name + "_InstAccParts_" + fn2
    
    print("Reading the File %s.dat" %(filename) )
    
    # Open the binary file and read the information within.
    p = open("%s/%s.dat" %(saved_inst_acc, filename), 'rb')
    IA_data = cP.load(p)
    p.close()
    
    # Return the information.
    return IA_data

def Read_InstAcc_Evolution(Cloud_name, resolution, tinit, tfin):
    
    global_data = {'info':'This dictionary contains the information of the cloud data as and its time evolution.'}

    global_data["time"] = []

    snapshot_data = Read_One_InstAcc(Cloud_name, resolution, 0)
    for key_name in snapshot_data.keys():
        global_data[key_name] = []

    # Number of snapshots.
    nsnaps = int((tfin - tinit)*10)

    for i in range(nsnaps+1):
        tnow = tinit + i/10.
        snapshot_data = Read_One_InstAcc(Cloud_name, resolution, i)
        global_data["time"].append(tnow)
        for key_name in snapshot_data.keys():
            if key_name != "tag":
                global_data[key_name].append(snapshot_data['%s'%key_name])
            
    return global_data



#################################################################################################
#  Read Global accretion.
#################################################################################################


def Read_One_GlobalAcc(Cloud_name, resolution, snapshot):
    """
    This function reads the binary where the cloud object information is stored.
    """

    import cPickle as cP

    # Set the string for the resolution.
    if resolution < 0.1:
        resolution_str = "%.3i" %(resolution*100)
    else:
        resolution_str = "%.2i" %(resolution*10)

    saved_inst_acc = "/data/gamera/jcibanezm/StratBox/AccretionPaper/KineticEnergy/GlobalAccretion/%s/%spc/" %(Cloud_name, resolution_str)

    snapshot_str      = str(snapshot)

    if (len(snapshot_str) == 1):
        snapshot_str = "0"+ snapshot_str
    fn2 = "snpsht"+snapshot_str

    # concatenate the filename.
    filename = Cloud_name + "_GlobalAccParts_" + fn2

    print("Reading the File %s.dat" %(filename) )

    # Open the binary file and read the information within.
    p = open("%s/%s.dat" %(saved_inst_acc, filename), 'rb')
    GA_data = cP.load(p)
    p.close()

    # Return the information.
    return GA_data

def Read_GlobalAcc_Evolution(Cloud_name, resolution, tinit, tfin):

    global_data = {'info':'This dictionary contains the information of the cloud data as and its time evolution.'}

    global_data["time"] = []

    snapshot_data0 = Read_One_GlobalAcc(Cloud_name, resolution, 0)
    for key_name in snapshot_data0.keys():
        global_data[key_name] = []

    # Number of snapshots.
    nsnaps = int((tfin - tinit)*10)

    for i in range(nsnaps+1):
        tnow = tinit + i/10.
        snapshot_data = Read_One_GlobalAcc(Cloud_name, resolution, i)
        global_data["time"].append(tnow)
        for key_name in snapshot_data0.keys():
            global_data[key_name].append(snapshot_data['%s'%key_name])

    return global_data

