def average_values_across_classes( unique_values_bonds, unique_values_angles, outputfilefolder ):
    #This function finds the average bond and angle term for each class. 
    
    ignore_rows = [] 

    #Find Average Values Bonds
    for i in range(0,(len(unique_values_bonds))):
        for j in range((i+1),len(unique_values_bonds)):
            #Finds if the bond class has already been encountered
            if ( unique_values_bonds[i][0] == unique_values_bonds[j][0] ) and (unique_values_bonds[i][1] ==  unique_values_bonds[j][1] )  or ( ( unique_values_bonds[i][0] == unique_values_bonds[j][1] ) and (unique_values_bonds[i][1] == unique_values_bonds[j][0] ) ):
                unique_values_bonds[i][2] = unique_values_bonds[i][2] + unique_values_bonds[j][2]
                unique_values_bonds[i][3] = unique_values_bonds[i][3] + unique_values_bonds[j][3]
                unique_values_bonds[i][4] = unique_values_bonds[i][4] + 1
                ignore_rows.append(j)
    
    #Average Bonds Printed
    fid = open(( outputfilefolder + 'Average_Modified_Seminario_Bonds' ), 'w')

    #Remove bond classes that were already present and find mean value
    for i in range(0,len(unique_values_bonds)):
        if (i in ignore_rows) == 0:
            unique_values_bonds[i][2] = unique_values_bonds[i][2] / unique_values_bonds[i][4]
            unique_values_bonds[i][3] = unique_values_bonds[i][3] / unique_values_bonds[i][4]
            fid.write( unique_values_bonds[i][0] + '-' + unique_values_bonds[i][1] + '  ' + str("%.2f" % unique_values_bonds[i][2]) + '  ' + str("%.3f" % unique_values_bonds[i][3]) + '\n' )

    fid.close()

    #Find Average Values Angles
    ignore_rows_angles = []
    
    #Find Average Values Angles
    for i in range(0,(len(unique_values_angles))):
        for j in range((i+1),len(unique_values_angles)):
            #Finds if the angle class has already been encountered
            if ( unique_values_angles[i][0] ==  unique_values_angles[j][0] and unique_values_angles[i][1] ==  unique_values_angles[j][1] and unique_values_angles[i][2] ==  unique_values_angles[j][2] ) or ( unique_values_angles[i][0] == unique_values_angles[j][2]  and unique_values_angles[i][1] == unique_values_angles[j][1] and unique_values_angles[i][2] == unique_values_angles[j][0] ):
                unique_values_angles[i][3] = unique_values_angles[i][3] + unique_values_angles[j][3]
                unique_values_angles[i][4] = unique_values_angles[i][4] + unique_values_angles[j][4]
                unique_values_angles[i][5] = unique_values_angles[i][5] + 1
                ignore_rows_angles.append(j)

    #Average Angles Printed
    fid = open(( outputfilefolder + 'Average_Modified_Seminario_Angles' ), 'w')

    #Remove angles classes that were already present and find mean value
    for i in range(0,len(unique_values_angles)):
        if (i in ignore_rows_angles ) == 0:
            unique_values_angles[i][3] = unique_values_angles[i][3] / unique_values_angles[i][5]
            unique_values_angles[i][4] = unique_values_angles[i][4] / unique_values_angles[i][5]
            fid.write( unique_values_angles[i][0] + '-' + unique_values_angles[i][1] + '-' + unique_values_angles[i][2] + '  ' + str("%.2f" % unique_values_angles[i][3]) + '  ' + str("%.3f" % unique_values_angles[i][4]) + '\n' )

    fid.close()
