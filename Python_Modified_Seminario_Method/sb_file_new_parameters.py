def sb_file_new_parameters(inputfilefolder, filename):
#Takes new angle and bond terms and puts them into a .sb file with name
#filename_seminario.sb

    fid_angles = open(inputfilefolder + 'Average_Modified_Seminario_Angles', 'r')
    fid_bonds = open(inputfilefolder + 'Average_Modified_Seminario_Bonds', 'r')

    angles = [] 
    bonds = []

    bond_lines = fid_bonds.readline()

    while bond_lines:
        bonds.append(bond_lines.strip().split('  '))
        bond_lines = fid_bonds.readline()

    angle_lines = fid_angles.readline()

    while angle_lines:
        angles.append(angle_lines.strip().split('  '))
        angle_lines = fid_angles.readline()

    #Opens Files
    fidout = open(inputfilefolder + filename  + '_Seminario.sb', 'wt') #Script produces this file
    fidout.write('*****                         Bond Stretching and Angle Bending Parameters - July 17*****\n')

    #Prints out bonds at top of file
    for i in range(0,len(bonds)):
        if (len(bonds[i]) == 5 ):
            fidout.write( bonds[i][0] + ' ' +  bonds[i][1] + '      ' + bonds[i][2] + '        Modified Seminario Method AEAA \n')
        else:
            fidout.write( bonds[i][0] + ' ' +  bonds[i][1] + '      ' +   bonds[i][2] + '        Modified Seminario Method AEAA \n' )

    fidout.write('\n')
    fidout.write('********                        line above must be blank\n')

    #Prints out angles in middle of file
    for i in range(0,len(angles)):
        if len(angles[i]) == 8:
            fidout.write( angles[i][0] + '    ' + angles[i][1] + '       ' +  angles[i][2] + '    Modified Seminario Method AEAA \n')
        else:
            fidout.write( angles[i][0] + '     ' +  angles[i][1] + '       ' +  angles[i][2] + '        Modified Seminario Method AEAA \n')

    fidout.write('\n')
    fidout.write('\n')

    fidout.close()
    fid_angles.close()
    fid_bonds.close()

