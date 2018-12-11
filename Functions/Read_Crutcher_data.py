
# coding: utf-8

# In[71]:

def Get_Crutcher_data():
    import numpy as np
    data_dir = "/home/jcibanezm/codes/StratBox/MagneticClouds"

    filename = "%s/apj333303t1_ascii.txt" %(data_dir)

    f = open(filename, "r")

    # Read the header
    Header = []
    for i in range(3):
        entry0 = f.readline()
        Header.append(entry0)

    # Read the names of the variables
    entry1 = f.readline()

    nH      = []
    Bz      = []
    sigma   = []
    Species = []
    ref     = []
    name    = []

    eof = False

    while not eof:   
    #for line in range(10):

        entry2 = f.readline().split("\t")

        if entry2[0] == "\n":
            eof = True
            break

        name.append(entry2[0])
        Species.append(entry2[1])
        ref.append(int(entry2[2]))

        exp     = int(entry2[3].split("^")[1])
        dens    =  float(entry2[3].split("x")[0])
        nh_here = dens*10**(exp)
        nH.append(nh_here)

        Bz.append(float(entry2[4]))
        sigma.append(float(entry2[5]))

    nH      = np.array(nH)
    Bz      = np.array(Bz)
    sigma   = np.array(sigma)
    Species = np.array(Species)
    ref     = np.array(ref)
    name    = np.array(name)

    f.close()

    Crutcher_data = {"info":"this dictionary contains the Zeeman data  information in Table 1, Crutcher et al 2010."}

    Crutcher_data["vars"] = entry1
    Crutcher_data["nH"] = nH
    Crutcher_data["Bz"] = Bz
    Crutcher_data["sigma"] = sigma
    Crutcher_data["ref"] = ref
    Crutcher_data["name"] = name
    Crutcher_data["Species"] = Species

    return Crutcher_data

