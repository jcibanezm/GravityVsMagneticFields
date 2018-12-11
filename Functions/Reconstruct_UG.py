
#######################################################

def Reconstruct_UG(main_cloud):
    """
    This function takes a cloud dictionary and reconstructs a Uniform Grid using the smallest resolution element in the dictionary as the UG resolution. 
    syntax:
        UG_object = Reconstruct_UG(cloud_dictionary)

    return:
        UG_object with all the same keys that where in the cloud_dictionary.
    
    """
    
    import numpy as np
    import cPickle as cP

    import Read_binary

    ppc = 3.0856776e18

    # Read the cloud data.
    #main_cloud = Read_binary.Read_One_Snapshot(Cloud_name, resolution, snapshot)

    # Make a UG of the main_cloud dictionary.
    # The Uniform grid is also a dictionary called Cloud_Grid.

    dxUG   = np.min(main_cloud["dx"])

    xminUG = np.min(main_cloud["x"]) - 4*dxUG
    xmaxUG = np.max(main_cloud["x"]) + 4*dxUG
    yminUG = np.min(main_cloud["y"]) - 4*dxUG
    ymaxUG = np.max(main_cloud["y"]) + 4*dxUG
    zminUG = np.min(main_cloud["z"]) - 4*dxUG
    zmaxUG = np.max(main_cloud["z"]) + 4*dxUG

    if xmaxUG < xminUG:
        print "The box goes around a periodic boundary. I'm transforming the x dimension by Delta x."
        positive = (main_cloud["x"]>=0)
        xx = main_cloud["x"]*positive + (xmaxUG + np.abs((xminUG - main_cloud["x"])))*np.logical_not(positive)
        
        main_cloud["x"] = xx
        xminUG = np.min(main_cloud["x"])
        xmaxUG = np.max(main_cloud["x"])


    LxUG = xmaxUG - xminUG 
    LyUG = ymaxUG - yminUG
    LzUG = zmaxUG - zminUG

    NxUG = int(LxUG / dxUG + 2)
    NyUG = int(LyUG / dxUG + 2)
    NzUG = int(LzUG / dxUG + 2)

    UG = np.zeros((NxUG, NyUG, NzUG))

    print ("===============================================")
    print ("Reconstructing a uniform grid for this cloud")
    print ("Uniform Grid properties:")
    print ("Lz = %.2f   Ly = %.2f   Lz = %.2f"%(LxUG/ppc, LyUG/ppc, LzUG/ppc))
    print ("dx = %.2f pc" %(dxUG/ppc))
    print ("Nx = %i,    Ny = %i,    Nz = %i" %(NxUG, NyUG, NzUG))
    print ("xmin = %.2f      xmax = %.2f" %(xminUG/ppc, xmaxUG/ppc))
    print ("ymin = %.2f      ymax = %.2f" %(yminUG/ppc, ymaxUG/ppc))
    print ("zmin = %.2f      zmax = %.2f" %(zminUG/ppc, zmaxUG/ppc))
    print ("===============================================")


    # Initialize a dictionary that will contain the grid where the cloud will be mapped.
    Cloud_Grid = {"hello": "I'm the new grid with the cloud information"}

    # Initialize basic properties of the uniform grid for the cloud. Reconstruct the x, y and z positions of the cloud uniform grid.
    # Also inizialize the cell size dx, dy, dz.

    X  = np.zeros_like(UG)
    DX = np.zeros_like(UG)

    for i in range(np.shape(UG)[0]):
        X[i,:,:]  = xminUG + i*dxUG
        DX[i,:,:] = dxUG

    Cloud_Grid["x"]  = X
    Cloud_Grid["dx"] = DX
    del X

    Y = np.zeros_like(UG)

    for j in range(np.shape(UG)[1]):
        Y[:,j,:]  = yminUG + j*dxUG
        DX[:,j,:] = dxUG
    
    Cloud_Grid["y"]  = Y
    Cloud_Grid["dy"] = DX
    del Y

    Z = np.zeros_like(UG)

    for k in range(np.shape(UG)[2]):
        Z[:,:,k]  = zminUG + k*dxUG
        DX[:,:,k] = dxUG
    
    Cloud_Grid["z"] = Z
    #Cloud_Grid["dz"] = DX
    del Z
    del DX

    ndens = np.zeros_like(UG)
    magx  = np.zeros_like(UG)
    magy  = np.zeros_like(UG)
    magz  = np.zeros_like(UG)
    velx  = np.zeros_like(UG)
    vely  = np.zeros_like(UG)
    velz  = np.zeros_like(UG)
    temp  = np.zeros_like(UG)
    gpot  = np.zeros_like(UG)

    for cc in range(len(main_cloud["x"])):
    
        DX2dx = int(main_cloud["dx"][cc]/dxUG)
    
        if DX2dx == 1:
            iindex   = int((main_cloud["x"][cc] - xminUG ) / dxUG + 0.1)
            jindex   = int((main_cloud["y"][cc] - yminUG ) / dxUG + 0.1 )
            kindex   = int((main_cloud["z"][cc] - zminUG ) / dxUG + 0.1 )
        
            ndens[iindex, jindex, kindex] = main_cloud["numdens"][cc]
        
        else:
            # I have to cycle Over the AMR grid.
            for sub_cycle in range(DX2dx**3):
                # Convert from a 1D array to a 3D array indexes.
                kk =  sub_cycle / (DX2dx * DX2dx)
                jj = (sub_cycle - kk*DX2dx*DX2dx) / DX2dx
                ii =  sub_cycle - kk*DX2dx*DX2dx - jj*DX2dx
            
                iindex   = int((main_cloud["x"][cc] - xminUG - main_cloud["dx"][cc]/2.0 + dxUG/2.) / dxUG + 0.1 ) + ii
                jindex   = int((main_cloud["y"][cc] - yminUG - main_cloud["dy"][cc]/2.0 + dxUG/2.) / dxUG + 0.1 ) + jj
                kindex   = int((main_cloud["z"][cc] - zminUG - main_cloud["dz"][cc]/2.0 + dxUG/2.) / dxUG + 0.1 ) + kk
            
                # Store the local density value here.
                ndens[iindex, jindex, kindex] = main_cloud["numdens"][cc]
                temp[iindex, jindex, kindex]  = main_cloud["temp"][cc]
            
                velx[iindex, jindex, kindex]  = main_cloud["velx"][cc]
                vely[iindex, jindex, kindex]  = main_cloud["vely"][cc]
                velz[iindex, jindex, kindex]  = main_cloud["velz"][cc]

                magx[iindex, jindex, kindex]  = main_cloud["magx"][cc]
                magy[iindex, jindex, kindex]  = main_cloud["magy"][cc]
                magz[iindex, jindex, kindex]  = main_cloud["magz"][cc]
                gpot[iindex, jindex, kindex]  = main_cloud["gpot"][cc]
    
    Cloud_Grid["numdens"] = ndens
    Cloud_Grid["temp"]    = temp

    Cloud_Grid["magx"] = magx
    Cloud_Grid["magy"] = magy
    Cloud_Grid["magz"] = magz
    Cloud_Grid["velx"] = velx
    Cloud_Grid["vely"] = vely
    Cloud_Grid["velz"] = velz
    Cloud_Grid["gpot"] = gpot

    del ndens 
    del magx  
    del magy  
    del magz  
    del velx  
    del vely  
    del velz
    del temp, gpot

    return Cloud_Grid
