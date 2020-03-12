# Assembler: 
import numpy as np

def load_mesh(in_file):
    vertices = np.genfromtxt(in_file+".ver").T
    with open (in_file+".ebv") as infile:
        inlist = infile.readlines()
    pre_list = [line.strip(' \n').split(',') for line in inlist]
    # CONTINUE HERE: 
    #        test if following line makes sense
    #        up to now, it is a list of arrays. Do we need an np.matrix filled up with 0s?
    elements = [np.array([int(st[k]) for k in range(len(st))]) for st in pre_list]

def assembler():
    # loop over elements and assemble stiffness matrix
    pass

def local_stiffness_tetra(nodes):
    """
    Input:  3x4 matrix, vertices of tetrahedron.
    Output: a vector with the coefficients of the local
            stiffness matrix. 
    """ 
    d = nodes.shape[0]
    stiff_rhs = np.vstack((np.zeros((1,d)),np.eye(d)))
    H = np.vstack((np.ones((1,d+1)),nodes))
    meas_t = np.absolute(np.linalg.det(H))/factorial(d)
    G = np.linalg.solve(H,stiff_rhs)
    stiffness_loc = meas_t*G@G.T
    return stiffness_loc.flatten()
        
def local_stiffness_prism(nodes):
    """    
    nodes: 3 x 6 matrix
    grad2D: matrix with rows (d alpha_i/dx, d alpha_i/dy)
    
    nodal functions are tensor products:

    basis: { phi_ij = alpha_i*beta_j}

    We identify which are the coordinates of the plane containing
    the triangles of the prism by subtracting adyacent nodes along "x3" axis
    and looking for the non-zero coordinate.
    
    A:  3 x 3 matrix of 2D integrals of the products between the
        2D gradientes of {alpha_i : i} with respect to the coordinates
        over the triangles of the prism.

    B:  2 x 2 matrix of 1D integrals of the products between the
        {beta_j : j} along the perpendicular coordinate. We use Simpson's
        rule.

    C:  3 x3 matrix of 2D integrals of the products between the
        {alpha_i : i}. We use Simpson's rule.

    D:  2 x 2 matrix of 1D integrals of the products between the
        1D derivatives of the {beta_j : j} along the perpendicular 
        coordinate.
    """
    difference = np.absolute(nodes[:,3]-nodes[:,0])
    axis = np.argmax(difference)
    triangle = np.setdiff1d(np.arange(3),[axis])
    stiff_rhs = np.vstack((np.zeros((1,2)),np.eye(2)))
    H1 = np.vstack((np.ones((1,3)),nodes[triangle,0:3]))
    grad2D = np.linalg.solve(H1,stiff_rhs)
    H2 = np.vstack((np.array([1,1]),nodes[axis,[0,3]]))
    deriv_axis = np.linalg.solve(H2,np.array([[0],[1]]))
    meas_t = np.linalg.norm(np.cross(nodes[:,1]-nodes[:,0],nodes[:,2]-nodes[:,0]))
    meas_l = difference[axis]
    A = meas_t*grad2D@grad2D.T
    B = meas_l*np.array([[1/3,1/6],[1/6,1/3]])
    C = meas_t/3*np.array([[1/2,1/4,1/4],[1/4,1/2,1/4],[1/4,1/4,1/2]])
    D = meas_l*deriv_axis@deriv_axis.T
    
    return np.kron(B,A) + np.kron(D,C)
    
def local_stiffness_pyramid(nodes):
    """
    nodes: 3 x 5 matrix
    centroid of a pyramid is a convex combination of base centroid and the top:
            (1/4)*Top + (1-1/4)*base_centroid


    TODO:   PONER POR QUE QUEDA submatriz identidad en G      

    """
    edges         = nodes[:,0,0,0,0,1,1,1,2,2,3]-nodes[:,1,2,3,4,2,3,4,3,4,4]
    diameter      = np.max(np.linalg.norm(edges,axis=0))
    base_centroid = np.mean(nodes[:,0:4],axis=1)
    centroid      = (nodes[:,4]+3*base_centroid)/4

    idea: 
    row_one       = np.hstack(([1], np.sum(nodes - centroid,axis=1)/(5*diameter)))
    lower_rows    = np.hstack((np.zeros(3,1),np.eye(3)))/diameter
    G             = np.vstack((row_one,lower_rows))

    return

local_stiffness = { 4 : local_stiffness_tetra,
                    5 : local_stiffness_pyramid,
                    6 : local_stiffness_prism}
