def dot_product(A, B):
    """
    Returns the dot product between two vectors.
    """

    import numpy as np

    A = np.array(A)
    B = np.array(B)

    AdotB = (A[0]*B[0] + A[1]*B[1] + A[2]*B[2])

    return AdotB



def dot_product_norm(A, B):
    """
    Returns the normalized dot product between two vectors.
    """

    import numpy as np

    A = np.array(A)
    B = np.array(B)
    
    # I should normalize vector B
    Bmag  = np.sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
    Amag  = np.sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2])    

    if (Bnorm == 0):
        B     = np.array([0., 0., 0.])
    else:
        B = B / Bmag    

    if Anorm == 0:
        A = np.array([0.,0.,0.])
    else:
        A = A / Amag

    AdotB = (A[0]*B[0] + A[1]*B[1] + A[2]*B[2])
    #AdotB = np.sqrt(AdotB*AdotB)
    
    return float(AdotB)


def cross_product_norm(A, B):
    """
    Returns the cross product between two vectors. Vector B is normalized.
    """

    import numpy as np

    A = np.array(A)
    B = np.array(B)
    
    # I should normalize vector B
    Bmag  = np.sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
    
    if np.min(Bmag) == 0:
        for i in range(len(Bmag)):
            offset = 1
            while Bmag[i] == 0:
                Bmag[i] = Bmag[i-offset]
                offset += 1
        
    Bnorm = B / Bmag
    B     = Bnorm
        
    Ci = A[1]*B[2] - A[2]*B[1]
    Cj = A[2]*B[0] - A[0]*B[2] 
    Ck = A[0]*B[1] - A[1]*B[0]
    
    C = np.array([Ci, Cj, Ck])
    Cmag = np.sqrt(Ci*Ci + Cj*Cj + Ck*Ck)
    
    return C, Cmag


def cross_product(A, B):
    """
    Returns the cross product between two vectors.
    """

    import numpy as np

    A = np.array(A)
    B = np.array(B)

    Ci = A[1]*B[2] - A[2]*B[1]
    Cj = A[2]*B[0] - A[0]*B[2]
    Ck = A[0]*B[1] - A[1]*B[0]

    C = np.array([Ci, Cj, Ck])
    Cmag = np.sqrt(Ci*Ci + Cj*Cj + Ck*Ck)

    return C, Cmag


##################################################

def Angle_vectors(vec1, vec2, rotate=False):
    """
    returns the relative angle between two vectors
    """
    import Vector_computations as vec
    import math
    import numpy as np
    
    # Compute the magnitude of vector 1.
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
    
    # Compute the dot product between the two unit vectors.
    vdotv_norm = vec.dot_product(v1norm, v2norm)
    
    # Calculate the angle between these two.
    angle_rad = np.arccos(vdotv_norm)
    
    # Convert these angles to degrees
    phi = angle_rad*180 / math.pi
    
    # Fold the angle distribution to be between 0-90 instead of 0-180
    for aa in range(len(phi)):
        delta_phi = 0
        if phi[aa] > 90:
            delta_phi = phi[aa] - 90
            phi[aa] = 90 - delta_phi
    
    # Rotate the relative angle distribution by 90 degrees.
    if rotate:
        phi = 90 - phi

    return phi 




def Relative_Angle_cos(vec1, vec2, rotate=False):
    """
    returns the relative angle between two vectors
    """
    import Vector_computations as vec
    import math
    import numpy as np

    # Compute the magnitude of vector 1.
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

    # Compute the dot product between the two unit vectors.
    vdotv_norm = vec.dot_product(v1norm, v2norm)

    return vdotv_norm

def Relative_Angle_atan(vec1, vec2, rotate=False):
    """
    returns the relative angle between two vectors
    """
    import Vector_computations as vec
    import math
    import numpy as np

    # Compute the magnitude of vector 1.
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

    
    # Get the dot product between these two vectors
    crossprod, crosnorm = cross_product(vec1, vec2)
    
    dotproduct = dot_product(vec1, vec2)
    # get the cross product between them

    phiangles = np.arctan(crosnorm/dotproduct)

    cosphi = np.cos(phiangles)
    
    return cosphi
