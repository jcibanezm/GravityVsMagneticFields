#!/home/jcibanezm/codes/ytAcc2/yt-x86_64/bin/python

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from yt.config import ytcfg;ytcfg["yt","__withinreason"]="True"

Cloud_name  = "M4e3"
resolution  = 0.06

offset      = 0
num_snaps   = 50

import yt
import os
import math
from yt import derived_field
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from yt.units import pc, kpc, second, Kelvin, gram, erg, cm
import cPickle as cP

import Vector_computations as vec

import timeit

# Define some constant parameters to be used.
mp      = 1.6726e-24  * gram # g
mu      = 1.2924
kb      = 1.3806e-16  *erg / Kelvin # erg K-1
GNewton = 6.6743e-8   * cm**3 / (gram * second**2 )# cm3 g-1 s-2
Msun    = 1.9884e33   * gram
mm      = mu*mp

ppc     = 3.0856776e18

# Create a derived field.
@derived_field(name="numdens", units="1/cm**3", force_override=True)
def numdens(field, data):
    dens_here = data["dens"].value
    dens_here = dens_here * gram / cm**3
    return dens_here/mm

yt.add_field('numdens', function=numdens, units="1/cm**3", force_override=True)

# Generalize for all the other clouds.
if Cloud_name == "M4e3":

    px = 180
    py = -25
    pz =  10
    cloud_alpha = 0.4

    rad = 40
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
    real_rad = 25

    px_str = str(px)
    py_str = str(py)
    pz_str = str(pz)

    px = px + 8
    py = py - 5
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


### HRO
def run_HRO(box, nmin, nmax, field):
    """
    Run the calculation of the HRO.
    """
    
    num_angles = 180

    angles = np.array(np.zeros(num_angles+1))
    hro    = np.array(np.zeros(num_angles))

    within_range = (box["numdens"]>=nmin) & (box["numdens"]<nmax) 

    vec1 = np.array([box["magnetic_field_x"]*within_range, box["magnetic_field_y"]*within_range, box["magnetic_field_z"]*within_range])
    #vec1 = np.array([box["density_gradient_x"]*within_range, box["density_gradient_y"]*within_range, box["density_gradient_z"]*within_range])

    if field == "velocity":
        
        BV   = box.quantities.bulk_velocity().in_cgs()
        velx = (box["velocity_x"].in_cgs() - BV[0])*within_range
        vely = (box["velocity_y"].in_cgs() - BV[1])*within_range
        velz = (box["velocity_z"].in_cgs() - BV[2])*within_range
        
        vec2 = np.array([velx, vely, velz])
        
    elif field == "density":
        
        vec2 = np.array([box["density_gradient_x"]*within_range, box["density_gradient_y"]*within_range, box["density_gradient_z"]*within_range])
    
    elif field == "gravity":
        
        vec2 = "do something"

    v1mag = np.sqrt(vec1[0]**2 + vec1[1]**2 + vec1[2]**2)

    # First check if any of the entries of vector 1 is equal to 0.
    # If true, delete those entries from both vector 1 and vector 2.    
    vec1_zeros = np.where(v1mag == 0)

    if len(vec1_zeros[0]) > 0:
        vec1x = np.delete(vec1[0], vec1_zeros)
        vec1y = np.delete(vec1[1], vec1_zeros)
        vec1z = np.delete(vec1[2], vec1_zeros)

        vec1  = np.array([vec1x, vec1y, vec1z])
        v1mag = np.sqrt(vec1[0]**2 + vec1[1]**2 + vec1[2]**2)

        vec2x = np.delete(vec2[0], vec1_zeros)
        vec2y = np.delete(vec2[1], vec1_zeros)
        vec2z = np.delete(vec2[2], vec1_zeros)

        vec2 = np.array([vec2x, vec2y, vec2z])
        v2mag = np.sqrt(vec2[0]**2 + vec2[1]**2 + vec2[2]**2)

    # Now check if from the remaining entries, there is any where vector 2 is zero.
    # if true, delete these entries as well.    
    v2mag = np.sqrt(vec2[0]**2 + vec2[1]**2 + vec2[2]**2)

    vec2_zeros = np.where(v2mag == 0)
    
    if len(vec2_zeros) > 0:    
        vec1x = np.delete(vec1[0], vec2_zeros)
        vec1y = np.delete(vec1[1], vec2_zeros)
        vec1z = np.delete(vec1[2], vec2_zeros)

        vec1  = np.array([vec1x, vec1y, vec1z])
        v1mag = np.sqrt(vec1[0]**2 + vec1[1]**2 + vec1[2]**2)


        vec2x = np.delete(vec2[0], vec2_zeros)
        vec2y = np.delete(vec2[1], vec2_zeros)
        vec2z = np.delete(vec2[2], vec2_zeros)

        vec2 = np.array([vec2x, vec2y, vec2z])
        v2mag = np.sqrt(vec2[0]**2 + vec2[1]**2 + vec2[2]**2)

    # Now, normalize the vectors
    v1norm = vec1 / v1mag
    v2norm = vec2 / v2mag

    # get the cross product between these two vectors
    crossprod, crosnorm = vec.cross_product(vec1, vec2)
    
    # Get the dot product between these two vectors
    dotproduct = vec.dot_product(vec1, vec2)
    
    # is the dot product normalized ? because I'm using crosnorm in the equation here!!!
    phiangles = np.arctan(crosnorm/dotproduct)

    # According to Juan, I need to take the cos of the angle here.
    # Because of the solid angle.
    cosphi = np.cos(phiangles)
    
    #cosphi = np.sin(phiangles)
    #cosphi = 180./(math.pi)*phiangles
    
    #print("min cosphi = ", np.min(cosphi), "max cosphi = ", np.max(cosphi))

    # range between 0,1 or -1, 1???
    dist = np.histogram(np.abs(cosphi), range=(0, 1), bins=num_angles, normed=True)

    hro     = dist[0]
    angles  = dist[1]
    
    # If I want to flip the array and get the angle with respect to the iso-density contour.
    if field == "density":
        hro = np.flipud(hro)
    
    return angles, hro

def calculate_xi(angles, hro):
    """
    Given the HRO array, compute the xi.
    
    Remember, angles are given in cos.phi.
    
    I should compute the sigma as well.
    """
    
    num_bins = np.shape(hro)[0]
    xi       = np.array(np.zeros(num_bins))
    sigma_xi = np.array(np.zeros(num_bins))
    
    for dindex in range(num_bins):
        # Perpendicular angles (given that the angles are the cos(phi). cos(90) = 0 and cos(0)=1. )
        for i in range(len(angles)):
            if angles[i] <= 0.25:
                index_cut_perp = i

        Aperp  = 0
        for i in range(index_cut_perp):
            area_here = 0
            area_here = (angles[i+1] - angles[i])*hro[dindex][i]
            Aperp  += area_here
            
        std_Aperp = np.std(hro[dindex][0:index_cut_perp])
        #print "$\\sigma_{\\parallel}$", std_Apar

        # Parallel angles
        for i in range(len(angles)):
            if angles[i] >= 0.75:
                index_cut_par = i
                break

        Apar  = 0
        for i in range(len(angles) - index_cut_par - 1):
            area_here = 0
            area_here = (angles[index_cut_par + i + 1] - angles[index_cut_par + i]) * hro[dindex][index_cut_par + i]
            Apar  += area_here
            
        std_Apar = np.std(hro[dindex][index_cut_par:len(angles)])
        #print "$\\sigma_{\\perp}$", std_Aperp

        xi_here = (Apar - Aperp) / (Apar + Aperp)
        
        sigma_xi2_here = 4.0 * (Aperp**2*std_Apar**2 + Apar**2*std_Aperp**2) / (Apar + Aperp)**4
        
        #print xi_here
        xi[dindex] = xi_here
        sigma_xi[dindex] = np.sqrt(sigma_xi2_here)
    
    return xi, sigma_xi


for snp in range(num_snaps):
    
    snapshot = snp + offset

    plt_file  = data_dir + basename + "plt_cnt_%.4i" %snapshot

    # Load data set
    pf = yt.load(plt_file)

    ### Slice of a cloud.
    # Move the center of the frame.
    px_now = px + bvx_mean * snapshot/10.
    py_now = py + bvy_mean * snapshot/10.
    pz_now = pz + bvz_mean * snapshot/10.

    cM  = [px_now*ppc, py_now*ppc, pz_now*ppc]
    rad = 80

    le = [cM[0] - rad*ppc, cM[1]-rad*ppc, cM[2]-rad/2.0*ppc]
    re = [cM[0] + rad*ppc, cM[1]+rad*ppc, cM[2]+rad/2.0*ppc]

    # Generate regions to access the particle and fluid data.
    box  = pf.region(cM, le, re)
    sph  = pf.sphere(cM, rad*pc)

    # why 23.0 ??
    slc = pf.slice("z", cM[2], center=cM, data_source=box)

    resolution_frb = 512
    #resolution_frb = 1024

    frb = slc.to_frb(rad*pc, resolution_frb, center=cM)["numdens"].value

    # Run the function in three different density ranges!!
    hros     = np.array([np.zeros(180), np.zeros(180), np.zeros(180)])
    hros_vel = np.array([np.zeros(180), np.zeros(180), np.zeros(180)])

    angles, hros[0] = run_HRO(box, 1.0e-2, 1.0,  "density")
    angles, hros[1] = run_HRO(box, 1.0,   1.0e2, "density")
    angles, hros[2] = run_HRO(box, 1.0e2, 1.0e4, "density")

    anglesv, hros_vel[0] = run_HRO(box, 1.0e-2, 1.0,  "velocity")
    anglesv, hros_vel[1] = run_HRO(box, 1.0,   1.0e2, "velocity")
    anglesv, hros_vel[2] = run_HRO(box, 1.0e2, 1.0e4, "velocity")

    # Logaritmic bins.
    numbins  = 20
    max_dens = math.log10(np.max(box["numdens"]))
    if max_dens > 4.30102999566: max_dens = 4.30102999566
    nbined   = np.logspace(-2, max_dens, num=numbins)

    long_hro     = [[0 for i in range(180)] for j in range(numbins)]
    long_hro_vel = [[0 for i in range(180)] for j in range(numbins)]

    for ii in range(numbins-1):
        # get nmin, nmax
        nmin, nmax = nbined[ii], nbined[ii+1]
        angles,  long_hro[ii]     = run_HRO(box, nmin, nmax,  "density")
        anglesv, long_hro_vel[ii] = run_HRO(box, nmin, nmax,  "velocity")

    xi, sigma   = calculate_xi(angles, long_hro)
    xiv, sigmav = calculate_xi(anglesv, long_hro_vel)

    xsize = 16  # 18
    ysize = 7.2 # 20

    vmin = 1.0e-3
    vmax = 1.0e4

    cmap = "RdYlBu_r"

    fig = plt.figure(figsize=(xsize, ysize))

    #ax0 = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False)
    ax0 = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=True)

    ax0.plot([0,1], [0,1], alpha=0)

    ax0.spines['right'].set_visible(False)
    ax0.spines['left'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)

    ax0.tick_params(axis="both", which="major", length=0, width=0, labelsize=0)

    #ax0.text(0.73, 0.13, "$\\xi = \, \\frac{A_{\\parallel} \,-\, A_{\\perp}}{A_{\\parallel} \,+\, A_{\\perp}} $" , fontsize=45)

    # ---------------------------------------------------------------------------------------
    # Cloud slice
    ax = fig.add_axes([0.039, 0.20, 0.32, 0.84])

    #extent = [cM[0]/ppc-rad, cM[0]/ppc+rad, cM[1]/ppc-rad, cM[1]/ppc+rad]
    extent = [-rad, rad, -rad, rad]

    cax = ax.imshow(frb, cmap=cmap, extent=extent, vmin=vmin, vmax=vmax, origin="lower",interpolation='gaussian', norm=LogNorm())

    # We change the fontsize of minor ticks label
    ax.tick_params(axis='both', which='major', length=10, width=2, labelsize=20)
    ax.tick_params(axis='both', which='minor', length=5, width=1.5, labelsize=12)


    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

    ax.set_xlabel("position x [pc]", fontsize=15, labelpad=0.0)
    ax.set_ylabel("position y [pc]", fontsize=15, labelpad=-2.0)

    ax.contour(frb, extent=extent, levels=[100], alpha=1, linewidths=2, colors="cyan")#colors='#7fc97f')
    ax.contour(frb, extent=extent, levels=[1],  alpha=1, linewidths=2,   colors="magenta")#colors='#beaed4')
    ax.contour(frb, extent=extent, levels=[0.01], alpha=1, linewidths=2, colors="k")#colors='#fdc086')

    # ------------------------------------------------------------------------------------
    # Colorbar
    cbar_ax = fig.add_axes([0.031, 0.12, 0.32, 0.075])
    cbar = fig.colorbar(cax, cax=cbar_ax, orientation='horizontal')

    cbar.set_label("number density [cm$^{-3}$]", fontsize=17)
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.tick_params('both', length=6, width=1.5, which='major')
    cbar.ax.tick_params('both', length=12, width=1.5, which='minor')
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.ax.xaxis.set_label_position("bottom")


    # ---------------------------------------------------------------------------------------
    # HRO - density gradient - magnetic field

    yymin = np.min(hros)
    yymax = np.max(hros)

    xxmin = np.min(np.arccos(angles)*180.0/math.pi)
    xxmax = np.max(np.arccos(angles)*180.0/math.pi)

    ax1 = fig.add_axes([0.42, 0.42, 0.25, 0.53])

    ax1.plot(np.flipud(np.arccos(angles[1:])*180.0/math.pi), np.flipud(hros[0]), "-", color="k",       linewidth=4)
    ax1.plot(np.flipud(np.arccos(angles[1:])*180.0/math.pi), np.flipud(hros[1]), "-", color="magenta", linewidth=4)
    ax1.plot(np.flipud(np.arccos(angles[1:])*180.0/math.pi), np.flipud(hros[2]), "-", color="cyan",    linewidth=4)

    ax1.set_xlabel("$\phi$", fontsize=20, labelpad=-3)
    ax1.set_ylabel("P(N)", fontsize=15)

    ax1.set_ylim(yymin*0.8, yymax*1.3)

    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

    ax1.tick_params(axis="both", which="major", length=8, width=2, labelsize=16)

    ax1.plot([22.5, 22.5], [-10, 10], "-k")
    ax1.plot([67.5, 67.5], [-10, 10], "-k")

    ax1.fill_between([0, 22.5],     [0, 0], [10, 10], color='#91bfdb', alpha='0.5')
    ax1.fill_between([67.5, 90.0],  [0, 0], [10, 10], color='#fc8d59', alpha='0.5')

    ax1.text(xxmax*0.1,  yymax*1.1,"$A_{\\parallel}$", fontsize=30)
    ax1.text(xxmax*0.8, yymax*1.1, "$A_{\\perp}$",    fontsize=30)

    ax1.text(xxmax*0.5, yymax*1.0,    "$10^{2}$ < n[cm$^{-3}$] <$10^{4}$",  fontsize=14, horizontalalignment='center', color="#00b2b2")
    ax1.text(xxmax*0.5, yymax*1.075,   "$10^{0}$ < n[cm$^{-3}$] <$10^{2}$",   fontsize=14, horizontalalignment='center', color="magenta")
    ax1.text(xxmax*0.5, yymax*1.15,     "$10^{-2}$< n[cm$^{-3}$] <$10^{0}$", fontsize=14, horizontalalignment='center', color="k")

    ax1.set_xlim(xxmin, xxmax)

    ax1.set_title("HRO $\mathrm{iso}\, n$-$\\hat{B}$", fontsize=20)

    # ---------------------------------------------------------------------------------------
    # xi plot
    ax2 = fig.add_axes([0.42, 0.08, 0.25, 0.25])

    ax2.set_xlabel("number density [cm$^{-3}$]", fontsize=15, labelpad=-2)
    ax2.set_ylabel("$\\xi_{nB}$", fontsize=17, labelpad=-5.0)

    ax2.plot([1.0e-2, 1.0e4], [0, 0], "-k", linewidth=1.0)

    ax2.plot(nbined, xi, "-ko", linewidth=2.0, markersize=5)
    ax2.errorbar(nbined, xi, yerr=sigma, fmt='o', color="k")

    ax2.fill_between([1.0e-2, 1.0e4], [0, 0], [1, 1],   color='#91bfdb', alpha='0.5')
    ax2.fill_between([1.0e-2, 1.0e4], [-1, -1], [0, 0], color='#fc8d59', alpha='0.5')

    ax2.set_ylim(-1.0, 1.0)
    ax2.set_xlim(1.0e-2, 1.0e4)

    ax2.set_xscale("log")

    # ---------------------------------------------------------------------------------------
    # HRO - velocity   - magnetic field

    yymin = np.min(hros_vel)
    yymax = np.max(hros_vel)

    ax3 = fig.add_axes([0.72, 0.42, 0.25, 0.53])

    ax3.plot(np.flipud(np.arccos(anglesv[1:])*180.0/math.pi), np.flipud(hros_vel[0]), "-", color="k",       linewidth=4)
    ax3.plot(np.flipud(np.arccos(anglesv[1:])*180.0/math.pi), np.flipud(hros_vel[1]), "-", color="magenta", linewidth=4)
    ax3.plot(np.flipud(np.arccos(anglesv[1:])*180.0/math.pi), np.flipud(hros_vel[2]), "-", color="cyan",    linewidth=4)

    ax3.set_xlabel("$\phi$", fontsize=20, labelpad=-3)
    ax3.set_ylabel("P(N)", fontsize=15)

    ax3.set_ylim(yymin*0.8, yymax*1.3)

    for tick in ax3.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax3.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

    ax3.tick_params(axis="both", which="major", length=8, width=2, labelsize=16)

    ax3.plot([22.5, 22.5], [-10, 10], "-k")
    ax3.plot([67.5, 67.5], [-10, 10], "-k")

    ax3.fill_between([0, 22.5],    [0, 0], [10, 10], color='#91bfdb', alpha='0.5')
    ax3.fill_between([67.5, 90.0], [0, 0], [10, 10], color='#fc8d59', alpha='0.5')

    ax3.text(8,  yymax*1.1, "$A_{\\parallel}$",        fontsize=30)
    ax3.text(72, yymax*1.1, "$A_{\\perp}$",    fontsize=30)

    ax3.text(45, yymax*1.0,   "$10^{2}$ < n[cm$^{-3}$] <$10^{4}$", fontsize=14, horizontalalignment='center', color="#00b2b2")
    ax3.text(45, yymax*1.075, "$10^{0}$ < n[cm$^{-3}$] <$10^{2}$", fontsize=14, horizontalalignment='center', color="magenta")
    ax3.text(45, yymax*1.15,  "$10^{-2}$< n[cm$^{-3}$] <$10^{0}$", fontsize=14, horizontalalignment='center', color="k")

    ax3.set_title("HRO $\\hat{v}$-$\\hat{B}$", fontsize=20)

    # ---------------------------------------------------------------------------------------
    # xi plot
    #ax2 = fig.add_axes([0.692, 0.43, 0.30, 0.42])
    ax4 = fig.add_axes([0.72, 0.08, 0.25, 0.25])

    ax4.set_xlabel("number density [cm$^{-3}$]", fontsize=15, labelpad=-2)
    ax4.set_ylabel("$\\xi_{vB}$", fontsize=17, labelpad=-5.0)

    ax4.plot([1.0e-2, 1.0e4], [0, 0], "-k", linewidth=1.0)

    ax4.plot(nbined,     xiv, "-ko", linewidth=2.0, markersize=5)
    ax4.errorbar(nbined, xiv, yerr=sigmav, fmt='o', color="k")

    ax4.fill_between([1.0e-2, 1.0e4], [0, 0], [1, 1],   color='#91bfdb', alpha='0.5')
    ax4.fill_between([1.0e-2, 1.0e4], [-1, -1], [0, 0], color='#fc8d59', alpha='0.5')

    ax4.set_ylim(-1.0, 1.0)
    ax4.set_xlim(1.0e-2, 1.0e4)

    ax4.set_xscale("log")

    for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(13)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(13)

    ax2.tick_params(axis="both", which="major", length=8, width=2)

    ax4 = fig.add_axes([0,0,1,1], frameon=False)

    ax4.text(0.001, 0.96, "a)", fontsize=25)
    ax4.text(0.375, 0.96, "b)", fontsize=25)
    ax4.text(0.69, 0.96,  "c)", fontsize=25)

    #fig.show()

    save_dir = "/data/gamera/jcibanezm/StratBox/MagneticCloudsPaper/Figure2"
    fig.savefig("%s/%s/%s_Slice_HRO_%.3i.pdf"%(save_dir, Cloud_name, Cloud_name, snapshot), format='pdf', dpi=100)
    print("Saving Figure %s/%s/%s_Slice_HRO_%.3i.pdf"%(save_dir, Cloud_name, Cloud_name, snapshot))
