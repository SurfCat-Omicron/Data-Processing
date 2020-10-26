import sys
import numpy
import glob
import os

"""Convert VAMAS file to TXT file.
Syntax:
    python convert_vms2txt.py [TYPE] [FILE]

TYPE: XPS or ISS - module used to import the VAMAS data (ISS default)
FILE: file to be converted. Will be renamed to the .txt ending."""
print('This script converts VAMAS to .txt (with header). Searches subdirectories and converts all found files.')
if len(sys.argv) > 2:
    index = 2
    if sys.argv[1].upper() == 'XPS':
        import XPS as module
        marker = 0
    elif sys.argv[1].upper() == 'ISS':
        import ISS as module
        marker = 1
else:
    index = 1
    marker = 1
    import ISS as module
    msg = 'No module specified. Using module ISS for import.'
    print(msg)  

getdir = os.getcwd()+'\\'+ input('\nDirectory:\n{}\\'.format(os.getcwd()))
os.chdir(getdir)
for filename in glob.glob('**/*.vms',recursive=True):
#filename = sys.argv[index]
    exp = module.Experiment(filename)

    for number in exp.cps.keys():
        print('Writing region {}'.format(number))
        new_filename = filename.split('--')[0] + '_Region_{}.txt'.format(number)
        f = open(new_filename, 'w+')
        f.writelines(['Energy,Counts/s\neV,A.U.\n'])
        if marker == 0: # xps
            x = exp.BE[number]
        elif marker == 1: # iss
            x = exp.energy[number]
        y = exp.cps[number]
        for i in range(len(y)):
            line = '{},{}\n'.format(x[i], y[i])
            f.write(line)
        f.close()
    print('Done!')
