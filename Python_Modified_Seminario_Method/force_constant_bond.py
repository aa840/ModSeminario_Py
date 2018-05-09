def force_constant_bond(atom_A, atom_B, eigenvalues, eigenvectors, coords):
    #Force Constant - Equation 10 of Seminario paper - gives force
    #constant for bond
    import numpy as np

    #Eigenvalues and eigenvectors calculated 
    eigenvalues_AB = eigenvalues[atom_A,atom_B,:]

    eigenvectors_AB = eigenvectors[:,:,atom_A,atom_B]

    #Vector along bond 
    diff_AB = np.array(coords[atom_B,:]) - np.array(coords[atom_A,:])
    norm_diff_AB = np.linalg.norm(diff_AB)

    unit_vectors_AB = diff_AB / norm_diff_AB
    
    k_AB = 0 

    #Projections of eigenvalues 
    for i in range(0,3):
        dot_product = abs(np.dot(unit_vectors_AB, eigenvectors_AB[:,i]))
        k_AB = k_AB + ( eigenvalues_AB[i] * dot_product )

    k_AB = -k_AB * 0.5 # Convert to OPLS form

    return k_AB






        
