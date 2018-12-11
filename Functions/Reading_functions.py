#!/home/jcibanezm/codes/yt26/yt-x86_64/bin/python

import os
import numpy as np

# Create an object with all the information.

class Cloud_properties(object):
    
    cloud_tag = ""
    cloud_mass, cloud_volume, cloud_cell_num, cloud_current_time                 = 0,0,0,0
    cloud_avg_dens, cloud_spherical_rad, cloud_tff                               = 0,0,0
    cloud_CM_x, cloud_CM_y, cloud_CM_z                                           = 0,0,0
    cloud_x_min, cloud_y_min, cloud_z_min, cloud_x_max, cloud_y_max, cloud_z_max = 0,0,0,0,0,0
    cloud_bulk_vel_x, cloud_bulk_vel_y, cloud_bulk_vel_z                         = 0,0,0
    cloud_x_vel_disp, cloud_y_vel_disp, cloud_z_vel_disp,  cloud_3D_vel_disp     = 0,0,0,0
    cloud_xy_area, cloud_yz_area, cloud_zx_area                                  = 0,0,0
    cloud_surf_dens_xy, cloud_surf_dens_yz, cloud_surf_dens_zx                   = 0,0,0
    cloud_avg_cs, cloud_sfr, cloud_sfr_JML, cloud_jeans_mass, cloud_jeans_yt     = 0,0,0,0,0
    cloud_particles, cloud_alpha_virial, cloud_EVT, cloud_kinetic_term           = 0,0,0,0
    cloud_gravitational_term, cloud_surface_pressure, cloud_magnetic_term        = 0,0,0    
    
    Header = ""
    
    # The class "constructor" - It's actually an initializer 
    def __init__(self, num_clouds, cloud_properties_names, data, Header = None):
        """
            Initialize the clouds properties object.
        """
        import re
        
        #for i in range(len(cloud_properties_names)):
        #    exec("self.%s = %s" %(cloud_properties_names[i], data[]))
        
        self.Header = Header
        
        for i in range(len(cloud_properties_names)):
            setattr(self, cloud_properties_names[i], data[i])


def read_cloud_parameters(sg_init, snapshot, filename="default", directory="MyCloudsProps", resolution=2, print_info=True):
    """
    input: 	sg_init    : time when self-gravity was turned on.
            snapshot     : number of snapshot.
            filename   : if filename does not follow the default naming, explicitly write the name of the file to read.
            directory  : Where are the files located.
            print_info : print some information to tell you how to distribute the cloud properties to individual arrays.
            
    Given the initial time of SG and the current snapshot, return the arrays containing the
    clouds physical parameters.
    Returns: 
            int   - num_clouds
            int   - num_parameters
            tuple - cloud physical parameters names
            array - Array of arrays containing the physical parameters
    """
    import numpy as np
    import os
    import math
    
    # Define some constant parameters to be used.
    kpc 	= 3.0856e21  # cm
    pc  	= 3.0856e18  # cm
    km  	= 1.0e5      # cm
    Myr 	= 3.1556e13  # s
    mp      = 1.6726e-24  # g
    mu      = 1.2924
    kb      = 1.3806e-16  # erg K-1
    GNewton = 6.6743e-8   # cm3 g-1 s-2
    Msun    = 1.9884e33   # g
    mm      = mu*mp
    
    # 
    fn0 = "CloudProps"
    fn1 = "t"+str(sg_init)
    current_snpshot = str(snapshot)
    if (len(current_snpshot) == 1):
        current_snpshot = "0"+ current_snpshot
    fn2 = "snpsht"+current_snpshot
    filename_here = fn0+"_"+fn1+"_"+fn2
    
    # name of the cloud properties file.
    # If the name is different than default, then give the filename explicitly.
    if (filename == "default"):
        name = filename_here+"_properties.dat"
    else:
        name = filename

    # Resolution is necessary for the directory name.
    if resolution > 1:
        res_str = str(int(resolution)) + "pc"
    else:
        res_str = "0" + str(int(resolution*10)) + "pc"
    
    # Directory where cloud properties are located.
    if (directory == "default"):
        directory_here = "MyCloudsProps/"+res_str
    else:
        directory_here = directory
    
    if not os.path.exists(directory_here):
        print "Not such directory, %s , are you sure the cloud properties are stored there ?" %(directory_here)
    
    # Basic variables.
    Header = []
    num_clouds = 0
    num_properties = 0
    
    # Open the file and start reading line by line.
    f = open("%s/%s"%(directory_here,name), "r")
    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())
    f.readline()                  # Blank space
    
    # Read the number of clouds. (I could put a conditional that I had indeed red the number of clouds)
    entry1     = f.readline().split("\t")
    num_clouds = int(entry1[1])
    
    # Read the number of properties.
    entry2     = f.readline().split("\t")
    num_properties = int(entry2[1])
    
    f.readline()                  # Blank space
    
    # Initialize the arrays of the names of the cloud properties and their values.
    properties_values      = []
    cloud_properties_names = ()
    
    # Read the names of the cloud properties and store it in this tuple.
    for i in range(num_properties):
        cloud_properties_names += (f.readline().split("\n")[0],)
    
    f.readline()                  # Blank space
    
    # Initialize an array for every entry in the cloud properties tuple.
    for prop_array in cloud_properties_names:
        exec("%s = []" % prop_array )
    
    # Loop over clouds and loop over properties.
    for i in range(num_clouds):
        for prop in cloud_properties_names:
            
            # There are two special properties that are not floating numbers, this properties are
            # the tag name and the number of cells in the cloud and are a string and an in respectively.
            # tore the value of the property in it's corresponding array.
            
            if prop == "cloud_tag":
                tag_here = [str(f.readline().split("\t")[1])]
                exec('{propname}.append({tag})'.format(propname=prop, tag='tag_here[0]'))
            elif prop == "cloud_cell_num":
                numcells = int(f.readline().split("\t")[1])
                exec("%s.append(%i)" % (prop, numcells))
            else:
                value_here = float(f.readline().split("\t")[1])
                exec("%s.append(%e)" % (prop, value_here))
    
    # Bookeeping: Make all arrays to be numpy arrays.
    for prop_array in cloud_properties_names:
        exec("%s = np.array(%s)" % (prop_array, prop_array))
    
    f.close()
    
    
    # Create a new file containing the properties + the snapshot number. So that I can load file after file
    # And each one will be saved in their own array name. The root of the name is given by the entries in
    # the tuple "cloud_properties_names" and the number of entries in each array is the number of clouds.
    
    snpshot_prop_names = ()
    for this_prop in cloud_properties_names:
        snpshot_prop_names += (this_prop+current_snpshot,)
        exec("%s = [] " % snpshot_prop_names[-1] )
    
    for i in range(num_clouds):
        j = 0
        for this_prop in cloud_properties_names:
            exec("value_here = %s[i]" %this_prop )
            if this_prop == "cloud_tag":
                value_here = value_here.split(" ")[0]
                exec('{propname}.append({tag})'.format(propname=snpshot_prop_names[j], tag='value_here'))
            elif this_prop == "cloud_cell_num":
                exec("%s.append(%s)" % (snpshot_prop_names[j], value_here))
            else:
                exec("%s.append(%e)" % (snpshot_prop_names[j], value_here))
            j +=1
    
    # Bookeeping: Make all arrays to be numpy arrays.
    array_of_arrays = []
    for this_prop in snpshot_prop_names:
        exec("%s = np.array(%s)" % (this_prop, this_prop))
        array_of_arrays.append(0)
        array_of_arrays[-1] = []
        for i in range(num_clouds):
            exec("array_of_arrays[-1].append(%s[i])" % this_prop )
    
    
    clouds_here = []
    clouds_here.append("num_clouds"+current_snpshot)
    clouds_here.append(num_clouds)
    
    props_here = []
    props_here.append("num_properties"+current_snpshot)
    props_here.append(num_properties)
    
    #return clouds_here, props_here, snpshot_prop_names, array_of_arrays
    return clouds_here, props_here, cloud_properties_names, array_of_arrays


#############################################################################
#
# 
#

def master_cloud_properties(sg_init_time, snapshot, filename="default", directory="default"):
    
    # Make a default filename given the sg_init_time, snapshot
    # Make a default folder name given the resolution, sg_init_time.
    
    a0, a1, a2, a3 = read_cloud_parameters(sg_init_time, snapshot, filename=filename, directory=directory, print_info=False)
    my_props = Cloud_properties(a0, a2, a3)
    
    return my_props


#############################################################################
#
# save the field information in the cloud object.
#

class Cloud_object(object):
    tag          = ""
    fields       = 0
    num_cells    = 0
    current_time = 0
    ncut         = 0
    Header = 0
    
    # The class "constructor" - It's actually an initializer 
    def __init__(self, tag, fields, num_cells, data, Header = None):
        """
            Initialize the cloud object. Set up the cloud tag, number of fields to be stored and 
            number of cells in the cloud.
        """
        import re
        import math
        
        # Define some constant parameters to be used.
        mp      = 1.6726e-24  # g
        mu      = 1.2924
        mm      = mu*mp
        
        self.tag = tag
        self.fields = fields
        self.num_cells = num_cells 
        self.Header = Header
        
        # Getting back the information from the cloud's tag.
        a0 = (tag.split('_')[0])
        a1 = (tag.split('_')[1])
        a2 = (tag.split('_')[2])
        
        tag0 = re.split('(\d+)', a0)
        tag1 = re.split('(\d+)', a1)
        tag2 = re.split('(\d+)', a2)
        
        self.current_time = tag1[1]
        self.ncut         = tag2[1]
        
        for i in range(len(fields)):
            setattr(self, fields[i], data[i])
        
        self.dens        = np.array(self.dens)
        self.Dens        = np.array(self.dens)
        self.Density     = np.array(self.dens)
        self.density     = np.array(self.dens)
        
        self.numdens     = np.array(self.dens) / mm
        self.ndens       = np.array(self.dens) / mm
        self.numDens     = np.array(self.dens) / mm
        self.numberDensity = np.array(self.dens) / mm
        self.numberdensity = np.array(self.dens) / mm
        
        self.temp        = np.array(self.temp)
        self.Temperature = np.array(self.temp)
        self.Temp        = np.array(self.temp)
        self.temperature = np.array(self.temp)
        
        self.velx        = np.array(self.velx)
        self.velocity_x  = np.array(self.velx)
        self.velocityx   = np.array(self.velx)
        
        self.vely        = np.array(self.vely)
        self.velocity_y  = np.array(self.vely)
        self.velocityy   = np.array(self.vely)
        
        self.velz        = np.array(self.velz)
        self.velocity_z  = np.array(self.velz)
        self.velocityz   = np.array(self.velz)
        
        velocity_magnitude = []
        for i in range(num_cells):
            vm2 = self.velx[i]*self.velx[i] + self.vely[i]*self.vely[i] + self.velz[i]*self.velz[i]
            vm = math.sqrt(vm2)
            velocity_magnitude.append(vm)
        
        self.velocity_magnitude = np.array(velocity_magnitude)
        self.vel_magnitude      = np.array(velocity_magnitude)
        self.vel_mag            = np.array(velocity_magnitude)

        self.magx        = np.array(self.magx)
        self.Bx          = np.array(self.magx)
        self.bx          = np.array(self.magx)
        self.B_x         = np.array(self.magx)
        self.b_x         = np.array(self.magx)
        self.magneticx   = np.array(self.magx)
        self.magnetic_x  = np.array(self.magx)
        
        self.magy        = np.array(self.magy)
        self.By          = np.array(self.magy)
        self.by          = np.array(self.magy)
        self.B_y         = np.array(self.magy)
        self.b_y         = np.array(self.magy)
        self.magneticy   = np.array(self.magy)
        self.magnetic_y  = np.array(self.magy)
        
        self.magz        = np.array(self.magz)
        self.Bz          = np.array(self.magz)
        self.bz          = np.array(self.magz)
        self.B_z         = np.array(self.magz)
        self.b_z         = np.array(self.magz)
        self.magneticz   = np.array(self.magz)
        self.magnetic_z  = np.array(self.magz)


#################################################################################
#
# given the name of the file, read the information inside and return an object
# containing the cloud fields.
#


def Get_cloud_info(cloud_filename, sg_init_time, snapshot, n_cut, resolution, dir_here='default'):
    """
    Read the cloud information stored in 'cloud_filename'. This will create an object containing the cloud's physical
    properties such as positions (x,y,z), cell size (dx, dy, dz), density, temperature, velocities (vx, vy, vz), 
    magnetic field strengths (magx, magy, magz), CellMass and CellVolume.
    
    cloud_filename = name of the file of the cloud I am interested in.
    sg_init_time   = time at which self-gravity was turned on (integer in Myr).
    snapshot       = snapshot that I want to read. Determines the evolutionary time.
    n_cut          = density threshold used to generate the cloud catalog.
    """
    # Read the cloud information.
    # Look for the directory where the cloud information was saved.
    
    print "---------------------------------------------------------------------------"
    print "Reading cloud's %s information ..." %cloud_filename
    
    
    # Directory where the information is stored
    if dir_here == 'default':
        directory_general = "./CloudObjects/"
    else:
        directory_general = dir_here
    
    if not os.path.exists(directory_general):
        print "The directory %s does not exist" %(directory_general)
    
    # Resolution of the data is part of the directory path.
    if resolution >= 1:
        resolution_str = str(int(resolution)) + "pc"
    else:
        resolution_str = "0"+str(int(resolution * 10)) + "pc"
    
    directory_resolution = directory_general+resolution_str+"/" 
    if not os.path.exists(directory_resolution):
        print "The directory %s does not exist" %(directory_resolution)
    
    # Time when self-gravity was turned on, snapshot and and density cut define the path as well.
    snapshot_str = str(snapshot)
    if len(snapshot_str) == 1 :
        snapshot_str = "0" + snapshot_str
    
    ncut_str = str(int(n_cut))
    
    directory_clouds = directory_resolution+"t"+str(sg_init_time)+"_snpshot"+snapshot_str+"_ncut"+ncut_str+"/"
    if not os.path.exists(directory_clouds):
        print "The directory %s does not exist" %(directory_clouds)
    
    # Final name of the file.
    temp_file = directory_clouds+cloud_filename
    
    
    # Read the data.
    Header = []
    # Open the file and start reading line by line.
    f = open(temp_file, "r")
    #f = open("%s/%s"%(directory_clouds,name_cloud), "r")
    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())
    Header.append(f.readline())
    f.readline()                  # Blank space
    
    # Read the number of clouds. (I could put a conditional that I had indeed red the number of clouds)
    entry0     = f.readline().split("\t")[1].split("\n")
    cloud_tag = entry0[0]
    
    # Read the number of clouds. (I could put a conditional that I had indeed red the number of clouds)
    entry1     = f.readline().split("\t")
    num_fields = int(entry1[1])
    
    # Read the number of properties.
    entry2     = f.readline().split("\t")
    num_cells = int(entry2[1])
    
    f.readline()                  # Blank space
    
    fields_line = f.readline().split("\t")
    fields = []
    for i in range(num_fields):
        fields.append(fields_line[i])
        
    data = []
    for i in range(num_fields):
        data.append(0)
        data[i] = []
        for j in range(num_cells):
            data[i].append(0)
    
    for j in range(num_cells):
        data_line = f.readline().split("\t")
        for i in range(num_fields):
            data[i][j] = float(data_line[i])
            
    f.close()
       
    data = np.array(data)
    
    # Save the cloud information in an object.
    my_cloud = Cloud_object(cloud_tag, fields, num_cells, data, Header=Header)
    
    print "Cloud %s stored in an object. The following fields are available:" %cloud_tag
    print fields
    
    return my_cloud

####################################################################
#
# returns a list with the names of the files containing the cloud 
# information
#

def master_read_cloud_information(sg_init_time=230, snapshot=0, n_cut=100, resolution=2, dir_here='default'):
    """
    this routine gives and array of objects containing the full cloud information, positions, densities, 
    temperatures, velocities, magnetic fiels, CellMass and CellVolume
    
    input parameters:
        sg_init_time : initial self-gravitating time
        snapshot     : snapshot at which I want to restore my information.
        ncut         : density threshold used to extract the cloud catalog.
        resolution   : resolution of the simulation used.
    
    """
    from os import walk
    
    # Directory where the information is stored
    if dir_here == 'default':
        directory_general = "./CloudObjects/"
    else:
        directory_general = dir_here

    if not os.path.exists(directory_general):
        print "The directory %s does not exist" %(directory_general)
    
    # Resolution of the data is part of the directory path.
    if resolution >= 1:
        resolution_str = str(int(resolution)) + "pc"
    else:
        resolution_str = "0"+str(int(resolution * 10)) + "pc"
    
    directory_resolution = directory_general+resolution_str+"/" 
    if not os.path.exists(directory_resolution):
        print "The directory %s does not exist" %(directory_resolution)
    
    # Time when self-gravity was turned on, snapshot and and density cut define the path as well.
    snapshot_str = str(snapshot)
    if len(snapshot_str) == 1 :
        snapshot_str = "0" + snapshot_str
    
    ncut_str = str(int(n_cut))
    
    directory_clouds = directory_resolution+"t"+str(sg_init_time)+"_snpshot"+snapshot_str+"_ncut"+ncut_str+"/"
    if not os.path.exists(directory_clouds):
        print "The directory %s does not exist" %(directory_clouds)
    
    cloud_filenames = []
    for (dirpath, dirnames, filenames) in walk(directory_clouds):
        cloud_filenames.extend(filenames)
        break
        
        
    cloud_list = []

    print "Number of cloud files here %i" %(len(cloud_filenames))

    for i in range(len(cloud_filenames)):
        cloud_list.append(Get_cloud_info(cloud_filenames[i], sg_init_time, snapshot, n_cut, resolution, dir_here))
        
    return cloud_list



###############################################################
#
# Object containing the re-computed cloud properties.
#

class regenerate_cloud_object(object):
        
    # The class "constructor" - It's actually an initializer 
    def __init__(self, tag, mass, volume, cell_num, avg_dens, spherical_rad, tff, CM_x, CM_y, CM_z,
                 x_min, x_max, y_min, y_max, z_min, z_max, bulk_vel_x, bulk_vel_y, bulk_vel_z,
                 bulk_vel_3D, x_vel_disp, y_vel_disp, z_vel_disp, vel_disp_3D, vel_disp_total, alpha_virial, avg_cs):
        """
            save the new cloud properties in an object.
        """
        
        self.cloud_tag      = tag
        self.cloud_mass     = mass
        self.cloud_volume   = volume
        self.cloud_cell_num = cell_num
        
        self.cloud_avg_dens = avg_dens
        self.cloud_spherical_rad = spherical_rad
        self.cloud_tff      = tff
        
        self.cloud_CM_x     = CM_x
        self.cloud_CM_y     = CM_y
        self.cloud_CM_z     = CM_z
        
        self.cloud_x_min    = x_min
        self.cloud_x_max    = x_max
        self.cloud_y_min    = y_min
        self.cloud_y_max    = y_max
        self.cloud_z_min    = z_min
        self.cloud_z_max    = z_max
        
        self.cloud_bulk_vel_x = bulk_vel_x
        self.cloud_bulk_vel_y = bulk_vel_y
        self.cloud_bulk_vel_z = bulk_vel_z
        
        self.cloud_bulk_vel_x = bulk_vel_x
        self.cloud_bulk_vel_y = bulk_vel_y
        self.cloud_bulk_vel_z = bulk_vel_z
        
        self.cloud_3D_vel_disp = vel_disp_3D
        self.cloud_3D_bulk_vel = bulk_vel_3D
        self.cloud_vel_disp_total = vel_disp_total
        self.cloud_avg_cs      = avg_cs
        self.cloud_alpha_virial= alpha_virial


#######################################################################################
#
# ...
#

def Regenerate_cloud_properties(sg_init_time, snapshot, resolution=2, n_cut=100, nmin=100, nmax=1.0e5, dir_here='default', object_here='cloud'):
    """
        get a list of clouds with the computed properties for the given density range.
        
        input parameters:
            sg_init_time : time when self-gravity was turned on.
            snapshot     : snapshot number.
            resolution   : resolution of the data used to generate the cloud catalog.
            n_cut        : density threshold used to generate the data.
            nmin         : minimum number density for the new estimation of the cloud parameters.
            nmax         : maximum number density for the new cloud properties.
	    dir_here     : path to the directory where the cloud objects are stored.
            object_here  : am I restoring the whole cloud or a clump inside the cloud. 'cloud' or 'clump'
    """
    
    import numpy as np
    import math
   
    print "Hey there I'm about to start doing stuff here"
    print "bla bla bla bla"
    print "..."
 
    # Define some constant parameters to be used.
    kpc = 3.0856e21  # cm
    pc  = 3.0856e18  # cm
    km  = 1.0e5      # cm
    Myr = 3.1556e13  # s
    
    mp      = 1.6726e-24  # g
    mu      = 1.2924
    kb      = 1.3806e-16  # erg K-1
    GNewton = 6.6743e-8   # cm3 g-1 s-2
    Msun    = 1.9884e33   # g
    mm      = mu*mp
     
    clouds_restored = master_read_cloud_information(sg_init_time=sg_init_time, snapshot = snapshot, n_cut=n_cut, resolution=resolution, dir_here=dir_here)
     
    cloud_list = []
    
    print "" 
    print "Hello there Juan, this is the number of clouds restored %i, at snapshot %i" %(len(clouds_restored), snapshot)
    
    for ii in range(len(clouds_restored)):
 
        ##################################################
        # Restore the cloud properties for a density range.
        ##################################################
                
        within_density_range = []
        above_density_threshold = []
        for j in range(len(clouds_restored[ii].CellMass)):
            if ( clouds_restored[ii].numdens[j] >= nmin) :
                above_density_threshold.append(1)
                if clouds_restored[ii].numdens[j] <= nmax :
                    within_density_range.append(1)
                else:
                    within_density_range.append(0)
            else:
                above_density_threshold.append(0)
                within_density_range.append(0)    
    
        within_density_range = np.array(within_density_range)
	above_density_threshold = np.array(above_density_threshold)        

        tag           = clouds_restored[ii].tag
        mass          = sum(clouds_restored[ii].CellMass * within_density_range)
        mass_total    = sum(clouds_restored[ii].CellMass )

        if mass ==  0:
            mass = 1.0e-99
        
        if object_here == 'clump':
            volume    = sum(clouds_restored[ii].CellVolume * above_density_threshold)
        elif object_here == 'cloud':
            volume    = sum(clouds_restored[ii].CellVolume)

        volume_total = sum(clouds_restored[ii].CellVolume)        
 
        cell_num      = sum(within_density_range) #len(clouds_restored[ii].dens)
        avg_dens      = mass_total / volume_total
        spherical_rad = (3 * volume / (4 * math.pi))**(1.0/3.0)
        tff           = math.sqrt( 3*math.pi / ( 32 * GNewton * avg_dens ) )
        CM_x          = sum(clouds_restored[ii].x * clouds_restored[ii].CellMass ) / mass
        CM_y          = sum(clouds_restored[ii].y * clouds_restored[ii].CellMass ) / mass
        CM_z          = sum(clouds_restored[ii].z * clouds_restored[ii].CellMass ) / mass
        
        x_min         = min(clouds_restored[ii].x * within_density_range)
        y_min         = min(clouds_restored[ii].y * within_density_range)
        z_min         = min(clouds_restored[ii].z * within_density_range)
        
        x_max         = max(clouds_restored[ii].x * within_density_range)
        y_max         = max(clouds_restored[ii].y * within_density_range)
        z_max         = max(clouds_restored[ii].z * within_density_range)
        
        bulk_vel_x    = sum(clouds_restored[ii].velx*clouds_restored[ii].CellMass * within_density_range) / mass
        bulk_vel_y    = sum(clouds_restored[ii].vely*clouds_restored[ii].CellMass * within_density_range) / mass
        bulk_vel_z    = sum(clouds_restored[ii].velz*clouds_restored[ii].CellMass * within_density_range) / mass
        bulk_vel_3D   = sum(clouds_restored[ii].velocity_magnitude*clouds_restored[ii].CellMass * within_density_range) / mass
        
        x_vel_disp    = math.sqrt(sum((clouds_restored[ii].velx - bulk_vel_x)**2 * clouds_restored[ii].CellMass * within_density_range) / mass )
        y_vel_disp    = math.sqrt(sum((clouds_restored[ii].velx - bulk_vel_x)**2 * clouds_restored[ii].CellMass * within_density_range) / mass )
        z_vel_disp    = math.sqrt(sum((clouds_restored[ii].velx - bulk_vel_x)**2 * clouds_restored[ii].CellMass * within_density_range) / mass )
        
        vel_disp_3D   = math.sqrt(sum((clouds_restored[ii].velocity_magnitude - bulk_vel_3D)**2 * 
                                      clouds_restored[ii].CellMass * within_density_range) / mass )
       
        avg_cs = sum( np.sqrt(clouds_restored[ii].temp *kb * 5.0 / 3.0 / mm ) * clouds_restored[ii].CellMass * within_density_range ) / mass
        
        vel_disp_total = math.sqrt(avg_cs**2 + 1.0/3.0 * vel_disp_3D**2)
         
	alpha_virial = 5 * vel_disp_total**2 * spherical_rad / ( GNewton * mass)

 
        # I need to compute the average sound speed of the cloud here.
        # Then I have to compute the total 1D velocity dispersion as sqrt(csound^2 + 1/3*velDisp_3D ).

        if mass > 0:

            print "***** Regenerating cloud object %i" %ii
            mycloud = regenerate_cloud_object(tag, mass, volume, cell_num, avg_dens, spherical_rad, tff, CM_x, CM_y, CM_z,
                                              x_min, x_max, y_min, y_max, z_min, z_max, bulk_vel_x, bulk_vel_y, bulk_vel_z,
                                              bulk_vel_3D, x_vel_disp, y_vel_disp, z_vel_disp, vel_disp_3D, vel_disp_total, 
      					      alpha_virial, avg_cs)
        
            cloud_list.append(mycloud)
    
    return cloud_list


def Regenerate_cloud_properties_new_notation(sg_init_time, snapshot, resolution=2, n_cut=100, nmin=100, nmax=1.0e5, dir_here='default', object_here='cloud'):
    """
        get a list of clouds with the computed properties for the given density range.

        input parameters:
            sg_init_time : time when self-gravity was turned on.
            snapshot     : snapshot number.
            resolution   : resolution of the data used to generate the cloud catalog.
            n_cut        : density threshold used to generate the data.
            nmin         : minimum number density for the new estimation of the cloud parameters.
            nmax         : maximum number density for the new cloud properties.
	    dir_here     : path to the directory where the cloud objects are stored.
            object_here  : am I restoring the whole cloud or a clump inside the cloud. 'cloud' or 'clump'
    """

    import numpy as np
    import math

    print "Hey there I'm about to start doing stuff here"
    print "bla bla bla bla"
    print "..."

    # Define some constant parameters to be used.
    kpc = 3.0856e21  # cm
    pc  = 3.0856e18  # cm
    km  = 1.0e5      # cm
    Myr = 3.1556e13  # s

    mp      = 1.6726e-24  # g
    mu      = 1.2924
    kb      = 1.3806e-16  # erg K-1
    GNewton = 6.6743e-8   # cm3 g-1 s-2
    Msun    = 1.9884e33   # g
    mm      = mu*mp

    clouds_restored = master_read_cloud_information(sg_init_time=sg_init_time, snapshot = snapshot, n_cut=n_cut, resolution=resolution, dir_here=dir_here)

    cloud_list = []

    print "" 
    print "Hello there Juan, this is the number of clouds restored %i, at snapshot %i" %(len(clouds_restored), snapshot)

    for ii in range(len(clouds_restored)):

        ##################################################
        # Restore the cloud properties for a density range.
        ##################################################

        within_density_range = []
        above_density_threshold = []
        for j in range(len(clouds_restored[ii].cell_mass)):
            if ( clouds_restored[ii].numdens[j] >= nmin) :
                above_density_threshold.append(1)
                if clouds_restored[ii].numdens[j] <= nmax :
                    within_density_range.append(1)
                else:
                    within_density_range.append(0)
            else:
                above_density_threshold.append(0)
                within_density_range.append(0)    

        within_density_range = np.array(within_density_range)
        above_density_threshold = np.array(above_density_threshold)        

        tag           = clouds_restored[ii].tag
        mass          = sum(clouds_restored[ii].cell_mass * within_density_range)
        mass_total    = sum(clouds_restored[ii].cell_mass )

        if mass ==  0:
            mass = 1.0e-99

        if object_here == 'clump':
            volume    = sum(clouds_restored[ii].cell_volume * above_density_threshold)
        elif object_here == 'cloud':
            volume    = sum(clouds_restored[ii].cell_volume)

        volume_total = sum(clouds_restored[ii].cell_volume)        

        cell_num      = sum(within_density_range) #len(clouds_restored[ii].dens)
        avg_dens      = mass_total / volume_total
        spherical_rad = (3 * volume / (4 * math.pi))**(1.0/3.0)
        tff           = math.sqrt( 3*math.pi / ( 32 * GNewton * avg_dens ) )
        CM_x          = sum(clouds_restored[ii].x * clouds_restored[ii].cell_mass ) / mass
        CM_y          = sum(clouds_restored[ii].y * clouds_restored[ii].cell_mass ) / mass
        CM_z          = sum(clouds_restored[ii].z * clouds_restored[ii].cell_mass ) / mass

        x_min         = min(clouds_restored[ii].x * within_density_range)
        y_min         = min(clouds_restored[ii].y * within_density_range)
        z_min         = min(clouds_restored[ii].z * within_density_range)

        x_max         = max(clouds_restored[ii].x * within_density_range)
        y_max         = max(clouds_restored[ii].y * within_density_range)
        z_max         = max(clouds_restored[ii].z * within_density_range)

        bulk_vel_x    = sum(clouds_restored[ii].velx*clouds_restored[ii].cell_mass * within_density_range) / mass
        bulk_vel_y    = sum(clouds_restored[ii].vely*clouds_restored[ii].cell_mass * within_density_range) / mass
        bulk_vel_z    = sum(clouds_restored[ii].velz*clouds_restored[ii].cell_mass * within_density_range) / mass
        bulk_vel_3D   = sum(clouds_restored[ii].velocity_magnitude*clouds_restored[ii].cell_mass * within_density_range) / mass

        x_vel_disp    = math.sqrt(sum((clouds_restored[ii].velx - bulk_vel_x)**2 * clouds_restored[ii].cell_mass * within_density_range) / mass )
        y_vel_disp    = math.sqrt(sum((clouds_restored[ii].velx - bulk_vel_x)**2 * clouds_restored[ii].cell_mass * within_density_range) / mass )
        z_vel_disp    = math.sqrt(sum((clouds_restored[ii].velx - bulk_vel_x)**2 * clouds_restored[ii].cell_mass * within_density_range) / mass )

        vel_disp_3D   = math.sqrt(sum((clouds_restored[ii].velocity_magnitude - bulk_vel_3D)**2 * 
                                      clouds_restored[ii].cell_mass * within_density_range) / mass )

        avg_cs = sum( np.sqrt(clouds_restored[ii].temp *kb * 5.0 / 3.0 / mm ) * clouds_restored[ii].cell_mass * within_density_range ) / mass

        vel_disp_total = math.sqrt(avg_cs**2 + 1.0/3.0 * vel_disp_3D**2)

        alpha_virial = 5 * vel_disp_total**2 * spherical_rad / ( GNewton * mass)


        # I need to compute the average sound speed of the cloud here.
        # Then I have to compute the total 1D velocity dispersion as sqrt(csound^2 + 1/3*velDisp_3D ).

        if mass > 0:

            print "***** Regenerating cloud object %i" %ii
            mycloud = regenerate_cloud_object(tag, mass, volume, cell_num, avg_dens, spherical_rad, tff, CM_x, CM_y, CM_z,
                                              x_min, x_max, y_min, y_max, z_min, z_max, bulk_vel_x, bulk_vel_y, bulk_vel_z,
                                              bulk_vel_3D, x_vel_disp, y_vel_disp, z_vel_disp, vel_disp_3D, vel_disp_total, 
      					      alpha_virial, avg_cs)

            cloud_list.append(mycloud)

    return cloud_list
