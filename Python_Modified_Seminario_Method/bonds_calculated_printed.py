def bonds_calculated_printed(outputfilefolder, vibrational_scaling_squared, bond_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords):
    #This function uses the Seminario method to find the bond
    #parameters and print them to file

    import numpy as np
    from force_constant_bond import force_constant_bond

    #Open output file bond parameters are written to 
    fid = open((outputfilefolder + 'Modified_Seminario_Bonds'), 'w')
    
    k_b = np.zeros(len(bond_list))
    bond_length_list = np.zeros(len(bond_list))
    unique_values_bonds = [] # Used to find average values 

    for i in range(0, len(bond_list)):
        AB = force_constant_bond(bond_list[i][0], bond_list[i][1],eigenvalues, eigenvectors, coords)
        BA = force_constant_bond(bond_list[i][1], bond_list[i][0],eigenvalues, eigenvectors, coords)
        
        # Order of bonds sometimes causes slight differences, find the mean
        k_b[i] = np.real(( AB + BA ) /2); 

        # Vibrational_scaling takes into account DFT deficities/ anharmocity    
        k_b[i] = k_b[i] * vibrational_scaling_squared

        bond_length_list[i] =  bond_lengths[bond_list[i][0]][bond_list[i][1]]
        fid.write(atom_names[bond_list[i][0]] + '-' + atom_names[bond_list[i][1]] + '  ')
        fid.write(str("%#.5g" % k_b[i])+ '   ' + str("%#.4g" % bond_length_list[i]) +  '   ' +
                  str(bond_list[i][0] + 1) +  '   ' + str(bond_list[i][1] + 1))
        fid.write('\n')

        unique_values_bonds.append([atom_names[bond_list[i][0]], atom_names[bond_list[i][1]], k_b[i], bond_length_list[i], 1 ])
    
    fid.close()

    return unique_values_bonds
