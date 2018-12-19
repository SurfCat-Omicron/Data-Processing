import sys
import numpy

"""Convert VAMAS file to CSV file.
Syntax:
    python convert_vms2csv.py [TYPE] [FILE]

TYPE: XPS or ISS - module used to import the VAMAS data (ISS default)
FILE: file to be converted. Will be renamed to the .csv ending."""

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
filename = sys.argv[index]
exp = module.Experiment(filename)

for number in exp.cps.viewkeys():
    print('Writing region {}'.format(number))
    new_filename = filename.split('--')[0] + '_Region_{}.csv'.format(number)
    f = open(new_filename, 'a')
    if marker == 0: # xps
        x = exp.BE[number]
    elif marker == 1: # iss
        x = exp.energy[number]
    y = exp.cps[number]
    for i in range(len(y)):
        line = '{},{}\r\n'.format(x[i], y[i])
        f.write(line)
    f.close()
print('Done!')
