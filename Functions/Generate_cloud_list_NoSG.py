#!/home/jcibanezm/codes/ytt30/yt-x86_64/bin/python

offset = 0

from yt.config import ytcfg;ytcfg["yt","__withinreason"]="True"
import os
import uuid

# Juan C. Ibanez-Mejia.   May 2015.   (AMNH/ITA)
# ---
# 
# Hi there, this is the notebook containing the algorithm used to extract the coud catalog of the stratified box simulations.
# This routine has a python executable counterpart which, is the one you should use to avoid having to run a python notebook to generate the cloud catalog, extract some properties of the clouds and also save the information of each cloud in an individual file.
# 
# This routine is intended to be run by manually changing the data file where the information is supposed to be extracted.
# It is intended to be machine independent, so depending on which machine I am running it (laptop or daikaiju) I select the appropriate path and files.

import yt
import numpy as np
import math
import pylab as P 
import os
from yt import derived_field

from yt.units import meter, gram, second, kilogram,  joule, cm, parsec, Kelvin, Megayear, kilometer, pc, kpc, km, erg

# Define some constant parameters to be used.
#kpc = 3.0856e21  # cm
#pc  = 3.0856e18  # cm
#km  = 1.0e5      # cm
#Myr = 3.1556e13  # s

ppc = 3.0856e18

mp      = 1.6726e-24  * gram # g
mu      = 1.2924
kb      = 1.3806e-16  * erg / Kelvin # erg K-1
GNewton = 6.6743e-8   * cm**3 / (gram * second**2 )# cm3 g-1 s-2
#Msun    = 1.9884e33   * gram 
mm      = mu*mp

machine = os.uname()[1]

# Create a derived field.
@derived_field(name="numdens", units="1/cm**3", force_override=True)
def numdens(field, data):
    #from yt.units import meter, gram, second, kilogram,  joule, cm, parsec, Kelvin, Megayear, kilometer, pc, kpc, km, erg
    #mp      = 1.6726e-24 * gram
    #mu      = 1.2924
    #mm      = mu*mp
    dens_here = data["dens"].value
    dens_here = dens_here * gram / cm**3
    return dens_here/mm

yt.add_field('numdens', function=numdens, units="1/cm**3", force_override=True )

################################################################################
######### Choose the proper path for the data depending on the machine i'm using
################################################################################
sg_init_time = 230
madir = "/data/manda/jcibanezm/StratBox/RealProd/2pc/1pc/NoSG_Evol/"

#snp_nr = 2411
#for qq in range(42):	
for qq in range(1):	
    
    snp_nr = (qq+offset)*10
    
    plt_file = "Strat_Box_1pc_NoSG_Particles_hdf5_plt_cnt_%.4i" %snp_nr
    par_file = "Strat_Box_1pc_NoSG_Particles_hdf5_part_%.4i"    %snp_nr
    	
    my_plt_file  = madir + plt_file
    my_part_file = madir + par_file
    
    # Load the data file
    pf = yt.load(my_plt_file, particle_filename=my_part_file)
    
    print ("=============================================================================")
    print ("Running the Get cloud list routine.")
    print ("Reading file: %s" %pf)
    print ("")
    
    
    # Define the cadence of the data files.
    current_time10 = int(pf.current_time.in_units("Myr").value * 10 + 2 )
    snapshot10= current_time10 - sg_init_time*10 + 8
    snapshot = snapshot10
    
    snapshot_str = str(snapshot)
    
    print("My Snapshot name is %s" %snapshot_str)
    
    # Define the region where I will extract the information.
    # Box +- 200 pc of the midplane
    c  = [  0.0*ppc, 0.0*ppc, 0.0*ppc]
    le = [ -500*ppc,-500*ppc,-200*ppc]
    re = [  500*ppc, 500*ppc, 200*ppc]
    
    box = pf.h.region(c,le,re)
    
    print("Generating object to extract clouds.")
    
    field = 'numdens'
    
    n_cut = 100 * cm**(-3)	
    c_max = box[field].max().in_cgs()
    
    # Do cloud cut to plot the cloud data later.
    print("Extracting connected contours.")
    cons, contours = box.extract_connected_sets(field, 1, n_cut, c_max)
    
    # Keep the clouds that are at least 100 cells in volume and drop the rest.
    num_contours = len(contours[0])
    fake_clouds  = 0
    all_clouds   = np.empty(num_contours)
    
    all_clouds.fill(1)
    
    
    j=0
    for i in range(num_contours-1):
        j = i+1
        obj = contours[0][j]["cell_mass"]
        if(obj.size < 100 ):
            fake_clouds   = fake_clouds + 1
            all_clouds[i] = -1
    
    # Real number of clouds & Create the final cloud_object.
    num_clouds = num_contours - fake_clouds
    cloud_obj = []
    
    # Fill the cloud_object array with each individual cloud.
    j=0
    for i in range(num_contours):
        if (all_clouds[i] != -1):
            #print "fake cloud", ii
            j = i + 1
            obj = contours[0][j]
            cloud_obj.append(obj)
    
    
    print (" Total number of contours found          : %i" %num_contours)
    print (" The number of tiny clouds is            : %i" %fake_clouds)
    print (" The total number of resolver clouds is  : %i" %num_clouds)
    
    # ### First Get the cloud properties.
    
    #################################################################
    # Set the name of the cloud folder with the snapshot time, 
    # time when self-gravity was turned on and the density threshold
    #################################################################
    
    # Get the filename to save the cloud information.
    # Filename: CloudsObject_t< time of init of SG >_snpsht < current snapshot >
    fn0 = "CloudsObject"
    fn1 = "t"+str(sg_init_time)
    
    if (len(snapshot_str) == 1):
        snapshot_str = "0"+ snapshot_str
    fn2 = "snpsht"+snapshot_str
    
    # concatenate the filename.
    filename = fn0+"_"+fn1+"_"+fn2
    
    # Define a tuple containing the various cloud physical properties, e.g. cloud_mass or cloud_velocity_dispersion
    # This tuple is going to initialize the arrays [].
    cloud_properties_names = "cloud_tag", "cloud_mass", "cloud_volume", "cloud_cell_num", "cloud_current_time", "cloud_avg_dens", \
    "cloud_spherical_rad", "cloud_tff", "cloud_CM_x",   "cloud_CM_y",     "cloud_CM_z", "cloud_x_min", "cloud_y_min", "cloud_z_min",\
    "cloud_x_max", "cloud_y_max", "cloud_z_max", "cloud_bulk_vel_x", "cloud_bulk_vel_y", "cloud_bulk_vel_z", "cloud_x_vel_disp", \
    "cloud_y_vel_disp", "cloud_z_vel_disp",  "cloud_3D_vel_disp", "cloud_xy_area", "cloud_yz_area", "cloud_zx_area",  "cloud_surf_dens_xy", \
    "cloud_surf_dens_yz", "cloud_surf_dens_zx", "cloud_avg_cs", "cloud_sfr", "cloud_sfr_JML", "cloud_jeans_mass", "cloud_num_particles", \
    "cloud_alpha_virial", "cloud_EVT", "cloud_kinetic_term", "cloud_gravitational_term", "cloud_surface_pressure", "cloud_magnetic_term"        
    
    # Initialize an aray [] for every entry in the tuple cloud_properties_names
    for prop_array in cloud_properties_names:
        exec("%s = []" % prop_array )
        
    particles_properties_names = "cloud_particles_tag", "cloud_particle_posx", "cloud_particle_posy", "cloud_particle_posz", \
    "cloud_particle_velx", "cloud_particle_vely", "cloud_particle_velz"
    
    for prop_array in particles_properties_names:
        exec("%s = []" % prop_array )
        
    
    print("Compute Cloud properties.")
    
    Debug = False
        
    #num_clouds = 2	
    for i in range(num_clouds):
        print("Running Cloud %i of %i" %(i, num_clouds))
             
        # ------------------------------------------------------------------------------------------------
        # Total Mass [Msun]
        # -----------------------------------------------------------------------------------------------
        tm = 0.0
        tM = np.sum(cloud_obj[i]["cell_mass"]).in_units("Msun")
        cloud_mass.append(tM)
    
        # ------------------------------------------------------------------------------------------------
        # Current time
        # ------------------------------------------------------------------------------------------------
        ct = 0
        ct = pf.current_time.in_units("Myr")
        cloud_current_time.append(ct)
        
        # ------------------------------------------------------------------------------------------------
        # Total Volume [cm3]
        # ------------------------------------------------------------------------------------------------
        tV = 0.0
        tV = np.sum(cloud_obj[i]["cell_volume"]).in_units("cm**3")
        cloud_volume.append(tV) 
        
        # ------------------------------------------------------------------------------------------------
        # Cloud cells [#]
        # ------------------------------------------------------------------------------------------------
        tCc = 0.0
        tCc = cloud_obj[i]["cell_mass"].size
        cloud_cell_num.append(tCc)    
        
        # ------------------------------------------------------------------------------------------------
        # Cloud average density  [g cm-3]
        # ------------------------------------------------------------------------------------------------
        avRho = 0.0
        avRho = cloud_mass[-1] / cloud_volume[-1]
        cloud_avg_dens.append(avRho)
    
        # ------------------------------------------------------------------------------------------------
        # Sperical radius. 
        # Assume the cloud's mass is concentrated in a spherical volume of uniform density [cm]
        # ------------------------------------------------------------------------------------------------
        Rad = 0.0
        Rad = ( cloud_volume[-1]*3.0 / (4.0 * math.pi) ) ** (1.0/3.0)
        cloud_spherical_rad.append(Rad)    
        
        # ------------------------------------------------------------------------------------------------
        # Cloud's free fall time. 
        # Assuming spherical cloud of uniform density [s]
        #          (   3 pi   )  1/2
        #    tff = (  -----   ) ^
        #          ( 32 G rho )
        # ------------------------------------------------------------------------------------------------
        tff = 1.0e99
        tff = math.sqrt( 3*math.pi / ( 32 * GNewton * cloud_avg_dens[-1] ) )
        cloud_tff.append(tff)
    
        # ------------------------------------------------------------------------------------------------
        # Cloud's center of mass [cm]
        # ------------------------------------------------------------------------------------------------
        cmx, cmy, cmz = cloud_obj[i].quantities.center_of_mass().in_units("cm")
        cloud_CM_x.append(cmx)
        cloud_CM_y.append(cmy)
        cloud_CM_z.append(cmz)
    
        # ------------------------------------------------------------------------------------------------
        # Cloud's bulk velocity [cm s-1]
        # ------------------------------------------------------------------------------------------------
        #bvx, bvy, bvz = cloud_obj[i].quantities.bulk_velocity().in_units("km/s")
        bvx, bvy, bvz = 0, 0, 0
        cloud_bulk_vel_x.append(bvx)
        cloud_bulk_vel_y.append(bvy)
        cloud_bulk_vel_z.append(bvz)
    
        # ------------------------------------------------------------------------------------------------
        # Cloud velocity dispersion [cm s-1]
        # ------------------------------------------------------------------------------------------------
        # z direction
        # Broken?
        #sigmaz = cloud_obj[i].quantities.weighted_variance("velz","cell_mass")[0]
        sigmaz = 0
        cloud_z_vel_disp.append(sigmaz)
    
        # y direciton
        #sigmay = cloud_obj[i].quantities.weighted_variance("vely","cell_mass")[0]
        sigmay = 0
        cloud_y_vel_disp.append(sigmay)
    
        # x dir
        #sigmax = cloud_obj[i].quantities.weighted_variance("velx","cell_mass")[0]
        sigmax = 0
        cloud_x_vel_disp.append(sigmax)
    
        # 3D velocity dispersion
        #sigma3D = cloud_obj[i].quantities.weighted_variance("velocity_magnitude","cell_mass")[0]
        sigma3D = 0
        cloud_3D_vel_disp.append(sigma3D)
    
        # ------------------------------------------------------------------------------------------------
        # Cloud's Jeans Mass smilarly computed as in Joung & Mac Low 2006 
        # Using a modified version of the sound speed assuming microturbulence contribution of the 
        # velocity dispersion.
        #
        #  Mj = rho_avg * lambdaJeans ^3
        #
        #  lambdaJeans = (pi / G rhoavg)^1/2 sigma_tot
        #  sigma_tot   = (cs^2 + 1/3 sigma3D^2 ) ^1/2
        # ------------------------------------------------------------------------------------------------
        # Cloud's average sound speed
        cs = 0.0
        #cs = cloud_obj[i].quantities.weighted_average_quantity("sound_speed","cell_mass")
        cloud_avg_cs.append(cs)
    
        sigma_tot   = np.sqrt( cloud_avg_cs[-1].value**2. + 1.0/3.0 * cloud_3D_vel_disp[-1].value**2. )
        lambdaJeans = np.sqrt( math.pi / (GNewton.value * cloud_avg_dens[-1].value) ) * sigma_tot
    
        # Jeans Mass
        Mjeans = cloud_avg_dens[-1] * lambdaJeans**3
        cloud_jeans_mass.append(Mjeans)
    
        # ------------------------------------------------------------------------------------------------
        # Cloud's alpha virial.  as in Bertoldi & McKee 1992
        #
        #            5 * sigma ^ 2 * R
        #   alpha =  ------------------
        #                 G * M
        #
        # Where sigma is the modified velocity dispersion (Chandrasekhar 1951) accounting for the thermal
        # and kinetic velocity dispersion (sigma^2 = (cs^2 + 1/3 sigma_kin^2))
        # ------------------------------------------------------------------------------------------------
        alphavir = 5 * sigma_tot **2 * cloud_spherical_rad[-1] / (GNewton * cloud_mass[-1] )    
        cloud_alpha_virial.append(alphavir)
        
        # ------------------------------------------------------------------------------------------------
        # Find the star formation  per free fall time of the cloud. As in Krumholz, Matzner & McKee 2006
        #                           -0.68        -0.32
        # SFR = 0.073 * alpha_vir ^       Mach ^ 
        # ------------------------------------------------------------------------------------------------
        Mach = cloud_3D_vel_disp[-1] / cloud_avg_cs[-1]
        sfr = 0.073 * cloud_alpha_virial[-1]**(-0.68) * Mach **(-0.32)
        Mdot = sfr * cloud_mass[-1] / (cloud_tff[-1] )
        
        # Star formation rate per free fall time of the cloud.
        cloud_sfr.append(Mdot)
        
        # Star formation rate as in Joung & Mac Low 2006.
        if (cloud_mass[-1] / cloud_jeans_mass[-1] > 1.0):
            sfr_JML = 0.3 * cloud_mass[-1] / ( cloud_tff[-1] )
        else:
            sfr_JML = 0.0
        cloud_sfr_JML.append(sfr_JML)
    
        # ------------------------------------------------------------------------------------------------
        # Tag the cloud with a unique name for each cloud. [string]
        # ------------------------------------------------------------------------------------------------
        # Cloud_number (tag0), Initial Cloud Time (tag1), DensityThresh (tag2)
        a0     = str(i)
        mydigs = len(a0)
        digs0  = 3
        if (mydigs < digs0):
            if (mydigs == 1):
                tag0 = "00" + a0
            if (mydigs == 2):
                tag0 = "0" + a0
        else:
            tag0 = a0
        
        tag1 = str(int((pf.current_time.in_units("Myr").value + 0.04 ) * 10 ))
        
        a2     = str(int(n_cut.value))
        mydigs = len(a2)
        digs2 = 3
        if (mydigs == 2 ):
            tag2 = "0" + a2
        else:
            tag2 = a2
                
        mytag = "c" + tag0 + "_t" + tag1 + "_n" + tag2
        cloud_tag.append(mytag)
            
        # ------------------------------------------------------------------------------------------------
        # Particles inside and particle tags [# , [strings(size #)]]
        # ------------------------------------------------------------------------------------------------
        # Put some if conditions asking if a particle file is present or not.
        c_sph  = [cmx.value, cmy.value, cmz.value]
        rad_sph = Rad
        
        if rad_sph <= box.index.get_smallest_dx():
            rad_sph = rad_sph*2
        
        sph = pf.sphere(c_sph, rad_sph)
        
        num_particles_inside = len(sph["particle_tag"])
            cloud_num_particles.append(num_particles_inside)
        
        cloud_particles_tag.append([])
        cloud_particle_posx.append([])
        cloud_particle_posy.append([])
        cloud_particle_posz.append([])
        cloud_particle_velx.append([])
        cloud_particle_vely.append([])
        cloud_particle_velz.append([])
        
        if cloud_num_particles > 0 :
            #cloud_particles_tag = []
            for prt_i in range(num_particles_inside):
                cloud_particles_tag[i].append(int(sph["particle_tag"][prt_i].value))
                cloud_particle_posx[i].append(sph["particle_posx"][prt_i].value)
                cloud_particle_posy[i].append(sph["particle_posy"][prt_i].value)
                cloud_particle_posz[i].append(sph["particle_posz"][prt_i].value)
                cloud_particle_velx[i].append(sph["particle_velx"][prt_i].value)
                cloud_particle_vely[i].append(sph["particle_velx"][prt_i].value)
                cloud_particle_velz[i].append(sph["particle_velx"][prt_i].value) 
           
        # ------------------------------------------------------------------------------------------------
        # NSurface Area and gas surface desity
        # ------------------------------------------------------------------------------------------------
        
        # First define a function to get the surface areas
        def get_surf_area(los, cmx, cmy, cmz, data_obj):
            field, weight_field = "dens", None
            weight_field = "d"+los
    
            if los=='x': perp1, perp2, losnum = 'dy', 'dz', 0
            if los=='y': perp1, perp2, losnum = 'dz', 'dx', 1
            if los=='z': perp1, perp2, losnum = 'dx', 'dy', 2
            
            # Project the cloud along a given line of sight
            prj = pf.proj(field, losnum, center=[cmx, cmy, cmz], weight_field=weight_field, data_source=data_obj)
    
            Sigma = 0*cm**2
            for i in range(len(prj["dens"])):
                if prj[los][i] == prj[los][i]:
                    sigmai = prj[perp1][i] * prj[perp2][i]
                    Sigma = Sigma + sigmai
            
            return Sigma
            
        #surf_xy = get_surf_area('z', cloud_CM_x[i], cloud_CM_y[i], cloud_CM_z[i], cloud_obj[i])
        #surf_yz = get_surf_area('x', cloud_CM_x[i], cloud_CM_y[i], cloud_CM_z[i], cloud_obj[i])
        #surf_zx = get_surf_area('y', cloud_CM_x[i], cloud_CM_y[i], cloud_CM_z[i], cloud_obj[i])
            
            surf_xy = 1
            surf_yz = 1
            surf_zx = 1	
    
        cloud_xy_area.append(surf_xy)
        cloud_yz_area.append(surf_yz)
        cloud_zx_area.append(surf_zx)
    
        # ------------------------------------------------------------------------------------------------
        # Cloud's mean Surface density [g cm-2]
        # ------------------------------------------------------------------------------------------------
        
        surf_dens_xy = cloud_mass[i] / cloud_xy_area[i]
        surf_dens_yz = cloud_mass[i] / cloud_yz_area[i]
        surf_dens_zx = cloud_mass[i] / cloud_zx_area[i]
    
        cloud_surf_dens_xy.append(surf_dens_xy)
        cloud_surf_dens_yz.append(surf_dens_yz)
        cloud_surf_dens_zx.append(surf_dens_zx)
    
        # ------------------------------------------------------------------------------------------------
        # Cloud's extention in x, y and z.
        # ------------------------------------------------------------------------------------------------
        xextrema = cloud_obj[i].quantities.extrema("x")
        yextrema = cloud_obj[i].quantities.extrema("y")
        zextrema = cloud_obj[i].quantities.extrema("z")
    
        cloud_x_min.append(xextrema[0].value)
        cloud_y_min.append(yextrema[0].value)
        cloud_z_min.append(zextrema[0].value)
        
        cloud_x_max.append(xextrema[1].value)
        cloud_y_max.append(yextrema[1].value)
        cloud_z_max.append(zextrema[1].value)
    
        # Need to work on this extensively !!!!
        # So far I have no computation of the Eulerian virial theorem parameters
    
        # ------------------------------------------------------------------------------------------------
        # Cloud's Eulerian Virial Theorem. As in McKee & Zweibel 1992
        #
        #    1 ..
        #    - I  = 2 (T - Ts) + Mg + W
        #    2      
        #                   1                                                 1
        # Kinetic E.   T =  - int  (3 Pth + rho vel^2) dV   surface P   Ts =  - int  Pth.r dS
        #                   2    V                                            2    S
        #                    1
        # Magnetic P   Mg =  --  int ( B^2 - B0^2 ) dV     "original McKee & Zweibel included the mag surface  
        #                   8 pi                            tension which I am ignoring here "
        #                  
        # Gravitational  W = int rho r.g dV
        #                   
        # ------------------------------------------------------------------------------------------------
        cloud_kinetic_term.append(0)
        cloud_gravitational_term.append(0)
        cloud_surface_pressure.append(0)
        cloud_magnetic_term.append(0)
        cloud_EVT.append(0)
    
        # Print some suff on screen.
        if Debug: 
            print "------------------------------------------------------------------"
            print "Calculated parameters of cloud   ", i
            print "Cloud tag                =",cloud_tag[-1]
            print "Mass                     =",cloud_mass[-1], "  [g]"
            print "Volume                   =",cloud_volume[-1] , "  [cm3]       ", cloud_volume[-1] / pc**3 , " [pc3]"
            print "Number of Cells          =",cloud_cell_num[-1]
            print "average density          =",cloud_avg_dens[-1], " [g/cm3]"
            print "Spherical radius         =",cloud_spherical_rad[-1], " [cm]     ", cloud_spherical_rad[-1] / pc, " [pc]" 
            print "Free fall time           =",cloud_tff[-1], " [s]     ", cloud_tff[-1] / Myr, " [Myr]" 
            print "Center of Mass           =",cloud_CM_x[-1],", ", cloud_CM_y[-1],", ", cloud_CM_z[-1], " [cm]"
    
            print "Bulk Velocity[cm/s]      =",cloud_bulk_vel_x[-1],", ", cloud_bulk_vel_y[-1],", ", cloud_bulk_vel_z[-1], " [cm/s]"
            print "Velocity disp[cm/s]      =",cloud_x_vel_disp[-1],", ", cloud_y_vel_disp[-1],", ", cloud_z_vel_disp[-1], " [cm s-1]"
            print "Velocity disp 3D         =",cloud_3D_vel_disp[-1], " [cm/s]"
    
            print "********************************** Jeans Mass Calculation: ***********************************"
            print "cs avg       =", cs ,         "     [cm/s]  "
            print "sigma turb   =", cloud_3D_vel_disp[-1], "     [cm/s]"
          
            print "sigma tot    =", sigma_tot,   "      [cm/s]    = ", sigma_tot / km, " [km/s] "
            print "jeans length =", lambdaJeans, "      [cm]      = ", lambdaJeans / pc," [pc] "
            print "Jeans Mass   =", cloud_jeans_mass[-1], "      [g] "
            
            print "Jeans Mass yt=", cloud_jeans_yt[-1],   "      Msun"
    
            print "********************************** Star Formation rate calculation ***********************************"
            print "alpha virial  =", cloud_alpha_virial[-1], "          Bertoldi & McKee 1992"
            print "Mach number   = ", Mach
            print "sfr_tff       = ", sfr                  , "          eq 41 in Krumholtz, Matzner & McKee 2006"
            print "Mdot          = ", cloud_sfr[-1]        , "          eq 42 ''     ''       ''        ''      "
            print ""
            print "SFR simple    = ", cloud_sfr_JML[-1], "          as in Joung & Mac Low 2006"
       
            print "********************************** Particles ***********************************"
            print "Only tracer particles present for now."
            print "Tracer Particles = ", cloud_particles[-1]
            print ""
            print "No calculation of:"
            print "                    the Euler Virial Theorem parameters."
            print "                    Cloud surface area." 
            print "                    Cloud surface density."
            print ""
            
    # Bookeeping: make clear that all this arrays are numpy arrays.
    for prop_array in cloud_properties_names:
        exec("%s = np.array(%s)" % (prop_array, prop_array))
    
    
    # ### Save a file with the cloud information.
    
    # I have to improve the format for the data that I am writing and also create an extra file for the particle tags for the different clouds.
    # I should actucally create an individual file containing the particle information in it.
    # 
    # particle position, particle position mapped in the grid, particle velocity, particle type (just in case the particle I'm looking at is not a tracer.)
    # 
    # And so on.
    
    # Make a directory where the cloud properties for each snapshot are to be stored.
    # Create a file with the initial SelfGravity time and the snapshot number in the name.
    # Save the various cloud properties in the file.
    
    dir_path = "/data/gamera/jcibanezm/StratBox/AccretionPaper/CloudPopulation/"
    
    # Get the filename to save the cloud information.
    # Filename: CloudsObject_t< time of init of SG >_snpsht < current snapshot >
    fn0 = "CloudProps"
    fn1 = "t"+str(sg_init_time)
    
    if (len(snapshot_str) == 1):
        snapshot_str = "0"+ snapshot_str
    fn2 = "snpsht"+snapshot_str
    
    # concatenate the filename.
    filename = fn0+"_"+fn1+"_"+fn2
    
    # Get the resolution for the folder name.
    small_dx = box.index.get_smallest_dx().in_units("pc").value
    if small_dx > 1:
        small_dx = str(int(small_dx+0.3))
    else:
        small_dx = "0" + str(int(small_dx*10+0.3))
    
    
    # Create path to the directory where info will be saved.
    directory = dir_path+"/MyCloudsProps/NoSG/"+small_dx+"pc"
    if not os.path.exists(directory):
        os.makedirs(directory)
    
        print("Saving Cloud properties: %s/%s" %(directory, filename+"_properties.dat"))
    
    name = filename+"_properties.dat"
    f = open("%s/%s"%(directory,name), 'w')
    print >> f, "# Header of the cloud list %s file." %filename
    print >> f, "# This file contains the list of clouds and some physical properties of the clouds."
    print >> f, "# Juan C. Ibanez-Mejia. 1st of March 2015.  @AMNH. \n"
    print >> f, "num_clouds\t", num_clouds
    print >> f, "num_properties\t", len(cloud_properties_names)
    print >> f, ""
    
    # Print the cloud name properties to allocate the arrays when I read the file.
    for prop in cloud_properties_names: 
        print >> f, prop, "\t",
    print >> f, ""
    
    # Write the name of the physical property (this_prop, e.g. cloud_mass) 
    # and it's corresponding value (value_here), separated by tabs.
    
    for i in range(num_clouds):
        for this_prop in cloud_properties_names:
            exec("value_here = %s[i]" % this_prop )
            print >> f, value_here,"\t",
        print >> f, ""
    f.close()
    
    
    # ### Save each cloud to an individual file.
    # Get the resolution for the folder name.
    small_dx = np.min(box['dx'].in_units("pc").value)
    if small_dx > 1:
        small_dx = str(int(small_dx+0.3))
    else:
        small_dx = "0" + str(int(small_dx*10+0.3))
    
    for ii in range(num_clouds):
    #for ii in range(1):
    
        # Create path to the directory where info will be saved.
        dir_path = "/data/gamera/jcibanezm/StratBox/AccretionPaper/CloudPopulation/"
        directory = dir_path + "/CloudObjects/NoSG/"+small_dx+"pc"
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        directory_snapshot = directory + "/t"+str(sg_init_time)+"_snpshot"+str(snapshot)+"_ncut"+str(int(n_cut.value))
        if not os.path.exists(directory_snapshot):
            os.makedirs(directory_snapshot)
        
        fields = ["x", "y", "z", "dx", "dy", "dz", "velx", "vely", "velz", "dens",  
                  "temp", "magx", "magy", "magz", "cell_mass", "cell_volume"]
            
        name_cloud = cloud_tag[ii]+".dat"
    
            print("Saving data of cloud %s/%s" %(directory_snapshot, name_cloud))
    
        f = open("%s/%s"%(directory_snapshot,name_cloud), 'w')
        print >> f, "# Header of the cloud %s object information." %cloud_tag[ii]
        print >> f, "# This file contains all the physical information of the cloud in case I need to retrieve it."
        print >> f, "# All data is in cgs units."
        print >> f, "# Juan C. Ibanez-Mejia. April 2015.  @AMNH & @ITA. \n"
        print >> f, "Cloud_name\t", cloud_tag[ii]
        print >> f, "num_fields\t", len(fields)
        print >> f, "num_cells\t", cloud_cell_num[ii]
        print >> f, ""
        
        ## Print the field names.
        for field_name in fields:
            f.write(field_name)
            f.write("\t")
        f.write("\n")
        
        # Print the information in the file.
        for j in range(cloud_cell_num[ii]):
            for k in range(len(fields)):
                print >> f, cloud_obj[ii][fields[k]][j].value, "\t",
            print >>f, "\n",
        f.close()
    
    
    # ### Save the tracer particles of each cloud to a file
    # Get the resolution for the folder name.
    small_dx = np.min(box['dx'].in_units("pc").value)
    if small_dx > 1:
        small_dx = str(int(small_dx+0.3))
    else:
        small_dx = "0" + str(int(small_dx*10+0.3))
    
    #for ii in range(num_clouds):
        for ii in range(0):
    
        # Create path to the directory where info will be saved.
        dir_path = "/data/gamera/jcibanezm/StratBox/AccretionPaper/CloudPopulation/"
        directory = dir_path + "/Particles_Clouds/NoSG/"+small_dx+"pc"
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        directory_snapshot = directory + "/t"+str(sg_init_time)+"_snpshot"+str(snapshot)+"_ncut"+str(int(n_cut.value))
        if not os.path.exists(directory_snapshot):
            os.makedirs(directory_snapshot)
        
        particles_fields = ["cloud_particles_tag", "cloud_particle_posx", "cloud_particle_posy", "cloud_particle_posz", 
                            "cloud_particle_velx", "cloud_particle_vely", "cloud_particle_velz" ]
            
        # When they become available I will use them.
        #"particle_cell_posx, "particle_cell_posy, "particle_cell_posz"
            
        name_cloud = cloud_tag[ii]+".dat"
        f = open("%s/%s"%(directory_snapshot,name_cloud), 'w')
        print >> f, "# Header of the cloud %s particles information." %cloud_tag[ii]
        print >> f, "# This file contains the particle tags, positions and velocities in a cloud."
        print >> f, "# All data is in cgs units."
        print >> f, "# Juan C. Ibanez-Mejia. Jan 2016.  @AMNH & @ITA. \n"
        print >> f, "Cloud_name\t", cloud_tag[ii]
        print >> f, "num_particle_properties\t", len(particles_fields)
        print >> f, "num_particles\t", cloud_num_particles[ii]
        print >> f, ""
        
        ## Print the field names.
        for field_name in particles_fields:
            print >> f, field_name, "\t",
        print >> f, ""
        
        for j in range(cloud_num_particles[ii]):
            for this_prop in particles_fields :
                exec("value_here = %s[ii][j]" % this_prop )
                print >> f, value_here,"\t",
            print >> f, ""
        
        f.close()
    
    
    # # This routine should end here !!!
    
