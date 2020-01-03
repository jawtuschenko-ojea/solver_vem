# Assembler: 
import numpy as np

local_stiffness = { 4 : local_stiffness_tetra,
                    5 : local_stiffness_pyramid,
                    6 : local_stiffness_prism}

def load_mesh(in_file):
    pass
    vertices = np.genfromtxt(in_file+".ver").T
    with open (in_file+."ebv") as infile:
        inlist = infile.readlines()
    pre_list = [line.strip(' \n').split(',') for line in inlist]
    # CONTINUE HERE: 
    #        test if following line makes sense
    #        up to now, it is a list of arrays. Do we need an np.matrix filled up with 0s?
    elements = [np.array([int(st[k]) for k in range(len(st))]) for st in pre_list]

def assembler():
    # loop over elements and assemble stiffness matrix
    pass

def local_stiffness_tetra(p):
    """
    Input:  3x4 matrix, vertices of tetrahedron.
    Output: a vector with the coefficients of the local
            stiffness matrix. 
    """ 
    d = p.shape[1]
    stiff_rhs = np.vstack((np.zeros((1,d)),np.eye(d)))
    H = np.vstack((np.ones((1,d+1)),nodes.T))
    meas_t = np.absolute(np.linalg.det(H))/factorial(d)
    G = np.linalg.solve(H,stiff_rhs)
    stiffness_loc = meas_t*G@G.T
    return stiffness_loc.flatten()
        
def local_stiffness_prism():
    pass

def local_stiffness_pyramid():
    pass