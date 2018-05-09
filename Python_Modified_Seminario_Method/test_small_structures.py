import os

directories = os.listdir('../Python_Testing_Final_Small_Structure/')

print(directories)

for dir in directories:
    print(dir)
    str = "python modified_Seminario_method.py '../Python_Testing_Final_Small_Structure/" + dir + "/' '../Python_Testing_Final_Small_Structure/" + dir + "/' 0.957"
    print(str)
    os.system(str)
