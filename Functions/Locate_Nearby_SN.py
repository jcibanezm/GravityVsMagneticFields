

#############################################
def read_SN_file():
    """
    This function reads the SN file and returns a dictionary with the time, positions and type of SN.
    """
    
    f = open('SNfeedback.dat', 'r')

    entry0 = f.readline()

    SN_dict = {"info":"this dictionary contains the information of the SN in the simulation."}

    num_props = len(entry0.split())
    props = []

    for i in range(num_props):
        props.append(str(entry0.split()[i].split("]")[1]))

    for i in range(num_props):
        SN_dict["%s" %(props[i])] = []

    eof = False

    while not eof:
        entry1 = f.readline()

        if entry1 == "":
            eof = True
            break

        for i in range(num_props):
            value_here = float(entry1.split()[i])
            SN_dict["%s" %(props[i])].append(value_here)

    f.close()
    
    return SN_dict


###############################################

def Locate_Nearby_SN(Cloud_name, min_dist=200):
    """
    Locate all the nearby SN explosions during the evolution of a given cloud.
    """
    
    import numpy as np
    
    if Cloud_name == "M4e3":
        px = 180
        py = -25
        pz =  7

        # Mean bulk velocity of the cloud.
        bvx_mean = -1.
        bvy_mean = -1.
        bvz_mean = -1.

    elif Cloud_name == "M3e3":

        px = 458
        py = -380
        pz = 17

        # Mean bulk velocity of the cloud.
        bvx_mean = 0.
        bvy_mean = 3.0
        bvz_mean = -2.0

    elif Cloud_name == "M8e3":
        px = 65
        py = 360
        pz = 25

        # Mean bulk velocity of the cloud.
        bvx_mean =  1.
        bvy_mean = -1.
        bvz_mean = -1.
    
    # Read the general SN file.
    SN_dict = read_SN_file()
    
    # Locate SN nearby.

    Nearby_SN = {"info":"this dictionary contains the SN that have exploded near the cloud"}

    Nearby_SN["time"] = []
    Nearby_SN["type"] = []
    Nearby_SN["distance"] = []

    Nearby_SN["px"] = []
    Nearby_SN["py"] = []
    Nearby_SN["pz"] = []


    ppc = 3.0856e18
    t0  = 7.54823618369e+15

    num_SN = len(SN_dict["time"])

    for ss in range(num_SN):

        t_now = SN_dict["time"][ss] - t0

        if t_now < 0:
            break

        # Get the position of the Cloud's center of mass.
        # units [cm]
        pxnow = px*ppc + bvx_mean*t_now
        pynow = py*ppc + bvy_mean*t_now
        pznow = pz*ppc + bvz_mean*t_now

        # distance to the SN.
        dx = SN_dict["posx"][ss] - pxnow
        dy = SN_dict["posy"][ss] - pynow
        dz = SN_dict["posz"][ss] - pznow

        dd = np.sqrt(dx**2 + dy**2 + dz**2)

        # I should account for the periodic boundary conditions.

        if dd < min_dist*ppc:

            #print SN_dict["radius"][ss]/ppc

            Nearby_SN["distance"].append(dd/ppc)
            Nearby_SN["time"].append(t_now / 3.1556e13)
            Nearby_SN["type"].append(SN_dict["type"][ss])

            Nearby_SN["px"].append(SN_dict["posx"][ss]/ppc - px)
            Nearby_SN["py"].append(SN_dict["posy"][ss]/ppc - py)
            Nearby_SN["pz"].append(SN_dict["posz"][ss]/ppc - pz)

    
    return Nearby_SN


### Print 

def Print_Nearby_SN(Nearby_SN):
    """
    Given the dictionary of the nearby SN explosions, print it nicely.
    """
    for i in range(len(Nearby_SN["time"])):
        print ("t = %.2f \t, d = %.2f \t, px = %.2f, py = %.2f, pz = %.2f" %(Nearby_SN["time"][i], Nearby_SN["distance"][i], \
               Nearby_SN["px"][i], Nearby_SN["py"][i], Nearby_SN["pz"][i]))
