#!/home/jcibanezm/codes/ytt30/yt-x86_64/bin/python
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from yt.config import ytcfg;ytcfg["yt","__withinreason"]="True"
import os
import uuid

Cloud_name  = "M3e3"
resolution  = 0.06
ncut        = 100
snapshot    = 20

# Use the offset if I have calculated the cloud properties for some snapshots already.
offset        = 0
num_snapshots = 50

import yt
import os
import math
from yt import derived_field
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from yt.units import pc, kpc, second, Kelvin, gram, erg, cm
import cPickle as cP

import timeit


# Define some constant parameters to be used.
mp      = 1.6726e-24  * gram # g
mu      = 1.2924
kb      = 1.3806e-16  *erg / Kelvin # erg K-1
GNewton = 6.6743e-8   * cm**3 / (gram * second**2 )# cm3 g-1 s-2
Msun    = 1.9884e33   * gram
mm      = mu*mp

ppc = 3.0856776e18

# Create a derived field.
@derived_field(name="numdens", units="1/cm**3", force_override=True)
def numdens(field, data):
    return data["dens"].in_cgs()/mm

yt.add_field('numdens', function=numdens, units="1/cm**3", force_override=True)

# Create all the derived fields calculating the force, per unit volume.


# Generalize for all the other clouds.
if Cloud_name == "M4e3":

    px = 180
    py = -25
    pz =  10
    cloud_alpha = 0.4

    rad = 50
    real_rad = 30

    px_str = str(px)
    py_str = str(py)
    pz_str = str(pz)

    pz = pz - 3.

    # Mean bulk velocity of the cloud.
    bvx_mean = -1.
    bvy_mean = -1.
    bvz_mean = -1.

    txt_xoffset = 15
    txt_yoffset = -25


if Cloud_name == "M3e3":

    px = 450
    py = -380
    pz = 25
    cloud_alpha = 0.4

    rad = 50
    real_rad = 30

    px_str = str(px)
    py_str = str(py)
    pz_str = str(pz)

    px = px + 8
    py = py
    pz = pz - 8

    # Mean bulk velocity of the cloud.
    bvx_mean = 0.
    bvy_mean = 3.0
    bvz_mean = -2.0

    txt_xoffset = 25
    txt_yoffset = -35

elif Cloud_name == "M8e3":
    px = 60
    py = 370
    pz = 30
    cloud_alpha = 0.8

    rad = 100
    real_rad = 30

    px_str = str(px)
    py_str = str(py)
    pz_str = str(pz)

    # Mean bulk velocity of the cloud.
    bvx_mean =  1.
    bvy_mean = -1.
    bvz_mean = -1.

    px = px + 5
    py = py - 10
    pz = pz - 5

    txt_xoffset = 30
    txt_yoffset = -50

cloud_dir = Cloud_name + "_" + "a%.2i" %(cloud_alpha*10) + "_x" + px_str + "_y" + py_str + "_z" + pz_str

if resolution < 0.1:
    resolution_str = "%.3i" %(resolution*100)
elif resolution == 1:
    resolution_str = "1"
else:
    resolution_str = "%.2i" %(resolution*10)


data_dir = "/data/manda/jcibanezm/StratBox/RealProd/1pc_and_AccClouds/AccPaper_Resims/Particles/" + cloud_dir + "/%spc/" %resolution_str

# Set the address and the name of the data set.
if resolution == 1:
    basename = "Strat_Box_" + Cloud_name + "_%spc_SG_hdf5_" %(resolution_str)
else:
    basename = "Strat_Box_" + Cloud_name + "_%spc_hdf5_" %(resolution_str)


for snp in range(num_snapshots):
    
    snapshot = snp + offset

    plt_file  = data_dir + basename + "plt_cnt_%.4i" %snapshot
    part_file = data_dir + basename + "part_%.4i"    %snapshot

    # Load data set
    pf = yt.load(plt_file , particle_filename=part_file)

    # Move the center of the frame.
    px_now = px + bvx_mean * snapshot/10.
    py_now = py + bvy_mean * snapshot/10.
    pz_now = pz + bvz_mean * snapshot/10.

    c = [px_now*ppc, py_now*ppc, pz_now*ppc]

    le = [c[0] - rad*ppc, c[1]-rad*ppc, c[2]-rad*ppc]
    re = [c[0] + rad*ppc, c[1]+rad*ppc, c[2]+rad*ppc]

    # Generate regions to access the particle and fluid data.
    box  = pf.region(c, le, re)
    sph  = pf.sphere(c, real_rad*pc)


    # In[3]:

    grad_fieldsx     = pf.add_gradient_fields(('flash', u'velx'))
    grad_fieldsy     = pf.add_gradient_fields(('flash', u'vely'))
    grad_fieldsz     = pf.add_gradient_fields(('flash', u'velz'))
    grad_fields_gpot = pf.add_gradient_fields(('flash', u'gpot'))


    # In[4]:

    def _conv_x(field, data):
        # Do I need to subtract the bulk velocity ?
        uxdx = data["velx"].in_cgs() * data["velx_gradient_x"].in_cgs()
        uxdy = data["vely"].in_cgs() * data["velx_gradient_y"].in_cgs()
        uxdz = data["velz"].in_cgs() * data["velx_gradient_z"].in_cgs()

        convx = uxdx + uxdy + uxdz
        return convx

    def _conv_y(field, data):
        uydx = data["velx"].in_cgs() * data["vely_gradient_x"].in_cgs()
        uydy = data["vely"].in_cgs() * data["vely_gradient_y"].in_cgs()
        uydz = data["velz"].in_cgs() * data["vely_gradient_z"].in_cgs()

        convy = uydx + uydy + uydz
        return convy

    def _conv_z(field, data):
        uzdx = data["velx"].in_cgs() * data["velz_gradient_x"].in_cgs()
        uzdy = data["vely"].in_cgs() * data["velz_gradient_y"].in_cgs()
        uzdz = data["velz"].in_cgs() * data["velz_gradient_z"].in_cgs()

        convz = uzdx + uzdy + uzdz
        return convz

    pf.add_field(('gas','conv_x'), function=_conv_x, units="cm/s**2", take_log=False,
                 display_name='convective term in the x direction.')

    pf.add_field(('gas','conv_y'), function=_conv_y, units="cm/s**2", take_log=False,
                 display_name='convective term in the y direction.')

    pf.add_field(('gas','conv_z'), function=_conv_z, units="cm/s**2", take_log=False,
                 display_name='convective term in the z direction.')


    # ---

    # In[5]:

    grad_fields_bx = pf.add_gradient_fields(('flash', u'magx'))
    grad_fields_by = pf.add_gradient_fields(('flash', u'magy'))
    grad_fields_bz = pf.add_gradient_fields(('flash', u'magz'))


    # In[6]:

    def _Lorentz_x(field, data):
        import math
        
        pi = math.pi
        term1 = (data["magx_gradient_z"] - data["magz_gradient_x"]) * data["magz"]
        term2 = (data["magy_gradient_x"] - data["magx_gradient_y"]) * data["magy"]
        
        force_x = 1.0 / (4.0 * pi) *(term1 - term2).in_cgs()
        return force_x

    def _Lorentz_y(field, data):
        import math
        
        pi = math.pi
        term1 = (data["magy_gradient_x"] - data["magx_gradient_y"]) * data["magx"]
        term2 = (data["magz_gradient_y"] - data["magy_gradient_z"]) * data["magz"]
        
        force_y = 1.0 / (4.0 * pi) *(term1 - term2).in_cgs()
        return force_y

    def _Lorentz_z(field, data):
        import math
        
        pi = math.pi
        term1 = (data["magz_gradient_y"] - data["magy_gradient_z"]) * data["magy"]
        term2 = (data["magx_gradient_z"] - data["magz_gradient_x"]) * data["magx"]
        
        force_z = 1.0 / (4.0 * pi) *(term1 - term2).in_cgs()
        return force_z


    # In[7]:

    pf.add_field(('gas','lorentz_x'), function=_Lorentz_x, units="g/(cm**2*s**2)", take_log=False, force_override=True,
                 display_name='Lorentz force in the x direction.', )

    pf.add_field(('gas','lorentz_y'), function=_Lorentz_y, units="g/(cm**2*s**2)", take_log=False, force_override=True,
                 display_name='Lorentz force in the y direction.', )

    pf.add_field(('gas','lorentz_z'), function=_Lorentz_z, units="g/(cm**2*s**2)", take_log=False, force_override=True,
                 display_name='Lorentz force in the z direction.', )


    # # Calculate the ram pressure and the rotation of the fluid.

    # In[8]:

    def _Vmag2(field, data):

        vmag2 = 0.5* (data["velx"] * data["velx"] + data["vely"] * data["vely"] + data["velz"] * data["velz"])
        return vmag2

    pf.add_field(('gas','vmag2'), function=_Vmag2, units="cm**2/s**2", take_log=False, force_override=True,
                 display_name='velocity magnitude squared.', )

    grad_fields_ram = pf.add_gradient_fields(('gas', 'vmag2'))


    # In[9]:

    grad_fields_vx = pf.add_gradient_fields(('flash', u'velx'))
    grad_fields_vy = pf.add_gradient_fields(('flash', u'vely'))
    grad_fields_vz = pf.add_gradient_fields(('flash', u'velz'))


    # In[10]:

    def _Vort_x(field, data):
        
        term1 = (data["velx_gradient_z"] - data["velz_gradient_x"]) * data["velz"]
        term2 = (data["vely_gradient_x"] - data["velx_gradient_y"]) * data["vely"]
        
        force_x = (term1 - term2).in_cgs()
        return force_x

    def _Vort_y(field, data):

        term1 = (data["vely_gradient_x"] - data["velx_gradient_y"]) * data["velx"]
        term2 = (data["velz_gradient_y"] - data["vely_gradient_z"]) * data["velz"]
        
        force_y = (term1 - term2).in_cgs()
        return force_y

    def _Vort_z(field, data):

        term1 = (data["velz_gradient_y"] - data["vely_gradient_z"]) * data["vely"]
        term2 = (data["velx_gradient_z"] - data["velz_gradient_x"]) * data["velx"]
        
        force_z = (term1 - term2).in_cgs()
        return force_z


    # In[11]:

    pf.add_field(('gas','curl_vort_x'), function=_Vort_x, units="cm/s**2", take_log=False, force_override=True,
                 display_name='Curl of the vorticity in the x direction.', )

    pf.add_field(('gas','curl_vort_y'), function=_Vort_y, units="cm/s**2", take_log=False, force_override=True,
                 display_name='Curl of the vorticity in the y direction.', )

    pf.add_field(('gas','curl_vort_z'), function=_Vort_z, units="cm/s**2", take_log=False, force_override=True,
                 display_name='curl of the vorticity in the z direction.', )


    # ---

    # In[12]:

    box  = pf.region(c, le, re)


    # In[13]:

    # Now, calculate here the advection on the reference frame of the box center of mass.


    # my idea is the following:
    # 
    # if I have a flow $V$, consisting on a background velocity $v_{0}$ and some fluctuations $\delta v$.
    # 
    # The advection term is, $$ (V \cdot \nabla)V = ((v_{0} + \delta v) \cdot \nabla) (v_{0} + \delta v) $$
    # 
    # given that $v_{0}$ is uniform, $\partial_{i} v_{0} = 0$, then $$(V \cdot \nabla)V  = ((v_{0} + \delta v) \cdot \nabla) \delta v $$ 
    # $$ (V \cdot \nabla)V = (v_{0} \cdot \nabla) \delta v  + (\delta v \cdot \nabla) \delta v $$
    # 
    # if I jump on an innertial reference frame, any alteration to the advective term will be contained in the $v_{0}$ term.
    # 
    # Therefore, in order to calculate the contribution of the advection, from the reference frame of the box containing my cloud, I should compute the following:
    # 
    # $$ (\delta v \cdot \nabla) \delta v =  (V \cdot \nabla)V - (v_{0} \cdot \nabla) \delta v $$
    # 
    # 

    # In[14]:

    def adv_newref_x(data):
        BV   = data.quantities.bulk_velocity().in_cgs()
        uxdx = (data["velx"].in_cgs() - BV[0] ) * data["velx_gradient_x"].in_cgs() 
        uxdy = (data["vely"].in_cgs() - BV[1] ) * data["velx_gradient_y"].in_cgs()
        uxdz = (data["velz"].in_cgs() - BV[2] ) * data["velx_gradient_z"].in_cgs()
        
        convx = uxdx + uxdy + uxdz
        return convx

    def adv_newref_y(data):
        BV   = data.quantities.bulk_velocity().in_cgs()
        uydx = (data["velx"].in_cgs() - BV[0] ) * data["vely_gradient_x"].in_cgs()
        uydy = (data["vely"].in_cgs() - BV[1] ) * data["vely_gradient_y"].in_cgs()
        uydz = (data["velz"].in_cgs() - BV[2] ) * data["vely_gradient_z"].in_cgs()

        convy = uydx + uydy + uydz
        return convy

    def adv_newref_z(data):
        BV   = data.quantities.bulk_velocity().in_cgs()
        uzdx = (data["velx"].in_cgs() - BV[0] ) * data["velz_gradient_x"].in_cgs()
        uzdy = (data["vely"].in_cgs() - BV[1] ) * data["velz_gradient_y"].in_cgs()
        uzdz = (data["velz"].in_cgs() - BV[2] ) * data["velz_gradient_z"].in_cgs()

        convz = uzdx + uzdy + uzdz
        return convz


    # ---

    # In[ ]:

    # Make a dictionary that contains the information.

    box_dict = {"info":"This dictionary contains the information in each cell of the box I am using to do the force balance and study the orientation of magnetic fields with several ISM properties."}

    box_dict = {"units":{"numdens":"1/cm**3", "cell_mass":"g", "dx":"cm", "temp":"K", "dens":"g/cm**3", "pres":"g/(cm*s**2)", 
                         "velx":"cm/s", "magx":"g**1/2/(cm**1/2*s)", "Pgrad_x":"g/(cm**2*s**2)",
                         "lorentz_x":"g/(cm**2*s**2)", "conv_x":"cm/s**2", "gacc_x":"cm/s**2", 
                         "ramp_x":"cm/s**2", "curl_vort_x":"cm/s**2"}}

    box_dict["numdens"]   = box["numdens"].in_cgs().value
    box_dict["cell_mass"] = box["cell_mass"].in_cgs().value
    box_dict["cell_volume"] = box["cell_volume"].in_cgs().value
    box_dict["dx"]        = box["dx"].in_cgs().value
    box_dict["temp"]      = box["temp"].in_cgs().value
    box_dict["dens"]      = box["dens"].in_cgs().value
    box_dict["pres"]      = box["pres"].in_cgs().value
    box_dict["velx"]      = box["velx"].in_cgs().value
    box_dict["vely"]      = box["vely"].in_cgs().value
    box_dict["velz"]      = box["velz"].in_cgs().value
    box_dict["magx"]      = box["magx"].in_cgs().value
    box_dict["magy"]      = box["magy"].in_cgs().value
    box_dict["magz"]      = box["magz"].in_cgs().value

    box_dict["Pgrad_x"]     = box["pressure_gradient_x"].in_cgs().value
    box_dict["Pgrad_y"]     = box["pressure_gradient_y"].in_cgs().value
    box_dict["Pgrad_z"]     = box["pressure_gradient_z"].in_cgs().value

    box_dict["lorentz_x"]   = box["lorentz_x"].in_cgs().value
    box_dict["lorentz_y"]   = box["lorentz_y"].in_cgs().value
    box_dict["lorentz_z"]   = box["lorentz_z"].in_cgs().value

    box_dict["conv_x"]      = box["conv_x"].in_cgs().value
    box_dict["conv_y"]      = box["conv_y"].in_cgs().value
    box_dict["conv_z"]      = box["conv_z"].in_cgs().value

    box_dict["conv_x2"]      = adv_newref_x(box).value
    box_dict["conv_y2"]      = adv_newref_y(box).value
    box_dict["conv_z2"]      = adv_newref_z(box).value

    box_dict["gacc_x"]      = box["gpot_gradient_x"].in_cgs().value
    box_dict["gacc_y"]      = box["gpot_gradient_y"].in_cgs().value
    box_dict["gacc_z"]      = box["gpot_gradient_z"].in_cgs().value

    box_dict["ramp_x"]      = box["vmag2_gradient_x"].in_cgs().value
    box_dict["ramp_y"]      = box["vmag2_gradient_y"].in_cgs().value
    box_dict["ramp_z"]      = box["vmag2_gradient_z"].in_cgs().value

    box_dict["curl_vort_x"] = box["curl_vort_x"].in_cgs().value
    box_dict["curl_vort_y"] = box["curl_vort_y"].in_cgs().value
    box_dict["curl_vort_z"] = box["curl_vort_z"].in_cgs().value


    # In[ ]:

    # Save the box information onto a pickle file.

    save_dir = "/data/gamera/jcibanezm/StratBox/MagneticCloudsPaper/PickleClouds/%s"%(Cloud_name)

    f = open('%s/Cloud_%s_snp%.3i_Forces.pkl'%(save_dir, Cloud_name, snapshot), 'wb') 
    cP.dump(box_dict, f, protocol=cP.HIGHEST_PROTOCOL)
    f.close()


# In[ ]:



