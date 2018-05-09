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

        #Gets Hessian 
        if len(tline) > 25 and (tline[0:24].strip() == 'Cartesian Force Constants'):
            tline = fid.readline()
            while len(tline) < 13 or (tline[0:12].strip() != 'Nuclear charges'):
                tmp = (tline.strip()).split()
                np.append(hessian, tmp, 0)
                tline = fid.readline()

            loop = 'n'

        tline = fid.readline()

    fid.close()

    list_coords = [float(x)*float(0.529) for x in list_coords]

    N = int( float(len(list_coords)) / 3.0 ) #Number of atoms

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

    return N
