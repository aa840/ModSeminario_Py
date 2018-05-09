def input_data_processing( inputfilefolder ):
    #This function takes input data that is need from the files supplied
    #Function extracts input coords and hessian from .fchk file, bond and angle
    #lists from .log file and atom names if a z-matrix is supplied

    import numpy as np
    import os.path


    #Gets Hessian in unprocessed form and writes .xyz file too 
    [unprocessed_Hessian, N, names, coords] =  coords_from_fchk( inputfilefolder, 'lig.fchk' );

    #Gets bond and angle lists
    [bond_list, angle_list] = bond_angle_list(inputfilefolder)

    with open("Number_to_Atom_type") as f:
        OPLS_number_to_name = f.readlines()

    OPLS_number_to_name = [x.split() for x in OPLS_number_to_name]; 

    length_hessian = 3 * N
    hessian = np.zeros((length_hessian, length_hessian))
    m = 0

    #Write the hessian in a 2D array format 
    for i in range (0,(length_hessian)):
        for j in range (0,(i + 1)):
            hessian[i][j] = unprocessed_Hessian[m]
            hessian[j][i] = unprocessed_Hessian[m]
            m = m + 1
  
    hessian = (hessian * (627.509391) )/ (0.529**2) ; #Change from Hartree/bohr to kcal/mol /ang

    # if zmat exists part here 
    atom_names = []

    for i in range(0,len(names)):
        atom_names.append(names[i].strip() + str(i + 1))

    if os.path.exists(inputfilefolder + 'Zmat.z'):    
        atom_names = []
        
        fid = open(inputfilefolder + 'Zmat.z') #Boss type Zmat
   
        tline = fid.readline()
        tline = fid.readline()

        #Find number of dummy atoms
        number_dummy = 0
        tmp = tline.split()
        
        while tmp[2] == '-1':
            number_dummy = number_dummy + 1
            tline = fid.readline()
            tmp = tline.split()

        if int(tmp[3]) < 800:
            for i in range(0,N):
                for j in range(0, len(OPLS_number_to_name)):
                    if OPLS_number_to_name[j][0] == tmp[3]:
                        atom_names.append(OPLS_number_to_name[j][1])
                
                tline = fid.readline()
                tmp = tline.split()
        else:
            #For CM1A format
            while len(tmp) < 2 or tmp[1] != 'Non-Bonded':
                tline = fid.readline()
                tmp = tline.split()
            
            tline = fid.readline()
            tline = fid.readline()
            tmp = tline.split()

            for i in range(0,N):
                atom_names.append(tmp[2])
                tline = fid.readline()
                tmp = tline.split()

        for i in range(0,N):
            if len(atom_names[i]) == 1:
                atom_names[i] = atom_names[i] + ' '
            
    return(bond_list, angle_list, coords, N, hessian, atom_names)

def coords_from_fchk(inputfilefolder, fchk_file):
    #Function extracts xyz file from the .fchk output file from Gaussian, this
    #provides the coordinates of the molecules
    import os.path
    import numpy as np

    if os.path.exists(inputfilefolder + fchk_file):
        fid = open((inputfilefolder + fchk_file), "r")
    else:
        fid_log = open((inputfilefolder + 'MSM_log'), "a")
        fid_log.write('ERROR = No .fchk file found.')
        fid_log.close()
        return 0,0
        
    tline = fid.readline()

    loop = 'y'
    numbers = [] #Atomic numbers for use in xyz file
    list_coords = [] #List of xyz coordinates
    hessian = []

    #Get atomic number and coordinates from fchk 
    while tline:
        #Atomic Numbers found
        if len(tline) > 16 and (tline[0:15].strip() == 'Atomic numbers'):
            tline = fid.readline()
            while len(tline) < 17 or (tline[0:16].strip() != 'Nuclear charges'):
                tmp = (tline.strip()).split()
                numbers.extend(tmp)
                tline = fid.readline()
            
        #Get coordinates
        if len(tline) > 31 and tline[0:31].strip() == 'Current cartesian coordinates':
            tline = fid.readline()
            while len(tline) < 15 or ( tline[0:14].strip() != 'Force Field' and tline[0:17].strip() != 'Int Atom Types'and tline[0:13].strip() != 'Atom Types'):
                tmp = (tline.strip()).split()
                list_coords.extend(tmp)
                tline = fid.readline()
            N = int( float(len(list_coords)) / 3.0 ) #Number of atoms

        #Gets Hessian 
        if len(tline) > 25 and (tline[0:26].strip() == 'Cartesian Force Constants'):
            tline = fid.readline()
            
            while len(tline) < 13 or (tline[0:14].strip() != 'Dipole Moment'):
                tmp = (tline.strip()).split()
                hessian.extend(tmp)
                tline = fid.readline()

            loop = 'n'

        tline = fid.readline()

    fid.close()

    list_coords = [float(x)*float(0.529) for x in list_coords]

    #Opens the new xyz file 
    file = open(inputfilefolder + 'input_coords.xyz', "w")
    file.write(str(N) + '\n \n')

    xyz = np.zeros((N,3))
    n = 0

    names = []

    fid_csv = open('elementlist.csv', "r")

    with fid_csv as f:
        lines = fid_csv.read().splitlines()

    element_names = []

    #Turn list in a matrix, with elements containing atomic number, symbol and name
    for x in range(0, len(lines)):
        element_names.append(lines[x].split(","))
        
    #Gives name for atomic number
    for x in range(0,len(numbers)):
        names.append(element_names[int(numbers[x]) - 1][1]) 

    #Print coordinates to new input_coords.xyz file
    for i in range(0, N):
        for j in range(0,3):
            xyz[i][j] = list_coords[n]
            n = n + 1

        file.write(names[i] + str(round(xyz[i][0],3)) + ' ' + str(round(xyz[i][1],3)) + ' ' + str(round(xyz[i][2], 3)) + '\n')

    file.close()
    return (hessian, N, names, xyz) 

def bond_angle_list(inputfilefolder):
#This function extracts a list of bond and angles from the Gaussian .log
#file
    import os.path

    fname = inputfilefolder + '/zmat.log'

    if os.path.isfile(fname) :
        fid = open(fname, "r")
    elif os.path.isfile(inputfilefolder + '/lig.log'):
        fid = open((inputfilefolder + '/lig.log'), "r")
    else:
        fid_log = open((inputfilefolder + 'MSM_log'), "a")
        fid_log.write('ERROR - No .log file found. \n')
        fid_log.close
        return 

    tline = fid.readline()
    bond_list = []
    angle_list = []

    n = 1
    n_bond = 1
    n_angle = 1
    tmp = 'R' #States if bond or angle
    B = []

    #Finds the bond and angles from the .log file
    while tline:
        tline = fid.readline()
        #Line starts at point when bond and angle list occurs
        if len(tline) > 80 and tline[0:81].strip() == '! Name  Definition              Value          Derivative Info.                !':
            tline = fid.readline()
            tline = fid.readline()
            #Stops when all bond and angles recorded 
            while ( ( tmp[0] == 'R' ) or (tmp[0] == 'A') ):
                line = tline.split()
                tmp = line[1]
                
                #Bond or angles listed as string
                list_terms = line[2][2:-1]

                #Bond List 
                if ( tmp[0] == 'R' ): 
                    x = list_terms.split(',')
                    #Subtraction due to python array indexing at 0
                    x = [(int(i) - 1 ) for i in x]
                    bond_list.append(x)

                #Angle List 
                if ( tmp[0] == 'A' ): 
                    x = list_terms.split(',')
                    #Subtraction due to python array indexing at 0
                    x = [(int(i) - 1 ) for i in x]
                    angle_list.append(x)

                tline = fid.readline()

            #Leave loop
            tline = -1

    return(bond_list, angle_list)

 
