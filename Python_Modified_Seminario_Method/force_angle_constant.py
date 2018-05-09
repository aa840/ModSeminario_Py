def  u_PA_from_angles( atom_A, atom_B, atom_C, coords ):
    #This gives the vector in the plane A,B,C and perpendicular to A to B

    import numpy as np
    from unit_vector_N import unit_vector_N

    diff_AB = coords[atom_B,:] - coords[atom_A,:]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB

    diff_CB = coords[atom_B,:] - coords[atom_C,:]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB

    u_N = unit_vector_N( u_CB, u_AB )

    u_PA = np.cross(u_N,  u_AB)
    norm_PA = np.linalg.norm(u_PA)
    u_PA = u_PA /norm_PA;

    return u_PA

def force_angle_constant( atom_A, atom_B, atom_C, bond_lengths, eigenvalues, eigenvectors, coords, scaling_1, scaling_2 ):
    #Force Constant- Equation 14 of seminario calculation paper - gives force
    #constant for angle (in kcal/mol/rad^2) and equilibrium angle in degrees

    import math
    import numpy as np 
    from unit_vector_N import unit_vector_N
    
    #Vectors along bonds calculated
    diff_AB = coords[atom_B,:] - coords[atom_A,:]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB

    diff_CB = coords[atom_B,:] - coords[atom_C,:]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB

    #Bond lengths and eigenvalues found
    bond_length_AB = bond_lengths[atom_A,atom_B]
    eigenvalues_AB = eigenvalues[atom_A, atom_B, :]
    eigenvectors_AB = eigenvectors[0:3, 0:3, atom_A, atom_B]

    bond_length_BC = bond_lengths[atom_B,atom_C]
    eigenvalues_CB = eigenvalues[atom_C, atom_B,  :]
    eigenvectors_CB = eigenvectors[0:3, 0:3, atom_C, atom_B]

    #Normal vector to angle plane found
    u_N = unit_vector_N( u_CB, u_AB )

    u_PA = np.cross(u_N,  u_AB)
    norm_u_PA = np.linalg.norm(u_PA)
    u_PA = u_PA / norm_u_PA

    u_PC = np.cross(u_CB, u_N)
    norm_u_PC = np.linalg.norm(u_PC)
    u_PC = u_PC / norm_u_PC

    sum_first = 0 
    sum_second = 0

    #Projections of eigenvalues
    for i in range(0,3):
        eig_AB_i = eigenvectors_AB[:,i]
        eig_BC_i = eigenvectors_CB[:,i]
        sum_first = sum_first + ( eigenvalues_AB[i] * abs(dot_product(u_PA , eig_AB_i) ) ) 
        sum_second = sum_second +  ( eigenvalues_CB[i] * abs(dot_product(u_PC, eig_BC_i) ) ) 

    #Scaling due to additional angles - Modified Seminario Part
    sum_first = sum_first/scaling_1 
    sum_second = sum_second/scaling_2

    #Added as two springs in series
    k_theta = ( 1 / ( (bond_length_AB**2) * sum_first) ) + ( 1 / ( (bond_length_BC**2) * sum_second) ) 
    k_theta = 1/k_theta

    k_theta = - k_theta #Change to OPLS form

    k_theta = abs(k_theta * 0.5) #Change to OPLS form

    #Equilibrium Angle
    theta_0 = math.degrees(math.acos(np.dot(u_AB, u_CB)))

    return k_theta, theta_0

def dot_product(u_PA , eig_AB):
    x = 0     
    for i in range(0,3):
        x = x + u_PA[i] * eig_AB[i].conjugate()
    return x 
