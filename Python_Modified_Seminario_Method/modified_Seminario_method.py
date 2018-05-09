def modified_Seminario_method(inputfilefolder, outputfilefolder, vibrational_scaling):
    #   Program to implement the Modified Seminario Method
    #   Written by Alice E. A. Allen, TCM, University of Cambridge
    #   Reference using AEA Allen, MC Payne, DJ Cole, J. Chem. Theory Comput. (2018), doi:10.1021/acs.jctc.7b00785

    # Pass as arguements the input folder containing the zmat.log/lig.log,
    # lig.fchk and optional Zmat.z file, the output folder where the new
    # parameters will be written and the vibrational frequency scaling constant
    # required.

    print(inputfilefolder)

    import time
    import numpy as np
    import os.path 

    from input_data_processing import input_data_processing
    from bonds_calculated_printed import bonds_calculated_printed
    from angles_calculated_printed import  angles_calculated_printed
    from average_values_across_classes import average_values_across_classes
    from sb_file_new_parameters import sb_file_new_parameters

    #Create log file 
    fid_log = open((inputfilefolder + 'MSM_log'), "w")
    fid_log.write('Modified Seminario Method \n')
    fid_log.write('Parametrization started for files in folder' + inputfilefolder + '\n')
    fid_log.write('Time is now: '+ time.strftime('%X %x %Z') + '\n')

    #Square the vibrational scaling used for frequencies
    vibrational_scaling_squared = vibrational_scaling**2; 

    #Import all input data
    [ bond_list, angle_list, coords, N, hessian, atom_names ] = input_data_processing( inputfilefolder)

    #Find bond lengths
    bond_lengths = np.zeros((N, N))

    for i in range (0,N):
        for j in range(0,N):
            diff_i_j = np.array(coords[i,:]) - np.array(coords[j,:])
            bond_lengths[i][j] =  np.linalg.norm(diff_i_j)

    eigenvectors = np.empty((3, 3, N, N), dtype=complex)
    eigenvalues = np.empty((N, N, 3), dtype=complex)
    partial_hessian = np.zeros((3, 3))

    for i in range(0,N):
        for j in range(0,N):
            partial_hessian = hessian[(i * 3):((i + 1)*3),(j * 3):((j + 1)*3)]
            [a, b] = np.linalg.eig(partial_hessian)
            eigenvalues[i,j,:] = (a)
            eigenvectors[:,:,i,j] = (b)

    # The bond values are calculated and written to file
    unique_values_bonds = bonds_calculated_printed( outputfilefolder, vibrational_scaling_squared, bond_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords )

    # The angle values are calculated and written to file
    unique_values_angles = angles_calculated_printed( outputfilefolder, vibrational_scaling_squared, angle_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords )

    #The final section finds the average bond and angle terms for each
    #bond/angle class if the .z exists to supply angle/bond classes and then 
    #writes the new terms to a .sb file
    if os.path.exists(inputfilefolder + 'Zmat.z'):
        average_values_across_classes( unique_values_bonds, unique_values_angles, outputfilefolder )
        sb_file_new_parameters(outputfilefolder, 'Python_Modified_Scaled')
        

import sys
modified_Seminario_method( sys.argv[1], sys.argv[2], float(sys.argv[3]) )
