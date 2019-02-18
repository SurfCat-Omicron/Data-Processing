import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

class Experiment(object):
    """Load an XPS experiment exported as text or VAMAS file.

Author: Jakob Ejler Sorensen
Version: 1.2
Date: 2017 July 13
    """

    def __init__(self, filename):
        """Initialize the class"""
        # Constants
        
        # Initialize variables
        self.KE = dict()
        self.BE = dict()
        self.cps = dict()
        self.dwell = dict()
        self.repeats = dict()
        self.mode = dict()
        self.mode_value = dict()
        self.line_energy = 0
        self.filekeys = {}
        self.first_comment = {}
        
        # Open filename with XPS data
        f = open(filename, 'r')
        print('Loading file: ' + filename + '\n')
        lines = f.readlines()
        f.close()
        
        # Read data from textfile:
        if filename.endswith('.txt'):
            raise ImportError('Text file import has not been implemented yet!')
        # Read data from VAMAS block file
        elif filename.endswith('.vms'):
            # Open filename with XPS data
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
	        # Old format:
            if lines[6].lower().startswith('experiment type'):
                self.format = 'Old VAMAS'
                # Find identifiers
                blocks_4 = [i for i in range(len(lines)) if (lines[i] == '-1\r\n') and (lines[i+1].lower() == 'kinetic energy\r\n')]
                blocks_2 = [i for i in range(len(lines)) if (lines[i] == 'XPS\r\n') and (lines[i+2] == '1253.6\r\n')]
                # Copy data points
                self.scans = len(blocks_4)
                for counter in range(self.scans):
                    i = blocks_4[counter]
                    self.mode[counter] = lines[i-11]
                    self.mode_value[counter] = float( lines[i-10] )
                    self.dwell[counter] = float( lines[i+9] )
                    self.repeats[counter] = float( lines[i+10] )
                    data_points = int( lines[i+16] )
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float( lines[i+4] )
                    E_start = float( lines[i+3] )
                    self.KE[counter] = np.arange(data_points)*E_step + E_start
                    self.BE[counter] = self.line_energy - self.KE[counter]
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = float( lines[i+19+counter_inner] )/self.dwell[counter]/self.repeats[counter]
            # New format
            elif lines[6].lower().startswith('created with'):
                self.format = 'New VAMAS'
                ending = filename.split('--')[1]
                ENDING = '--{}_{}-Detector_Region.vms'
                filen = filename.rstrip(ending)
                self.filename = filen

                # Detect number of regions in XPS data
                for counter in range(100):
                    for subcounter in range(1000):
                        try:
                            f = open(filen + ENDING.format(counter+1, subcounter+1), 'r')
                            if not counter in self.filekeys:
                                self.filekeys[counter] = []
                            self.filekeys[counter].append(subcounter)
                            f.close()
                        except IOError:
                            break
                    if subcounter == 0:
                        break

                # Open filename with XPS data
                for counter in self.filekeys.keys():
                    for subcounter in self.filekeys[counter]:
                        new_filename = filen + ENDING.format(counter+1,subcounter+1)
                        f = open(new_filename, 'r')
                        lines = f.readlines()
                        f.close()
                        print('Loading file: ' + new_filename)
                        print('New format')
                        # ID block 4 by their "-1" line
                        blocks_4 = [i for i in range(len(lines)) if (lines[i] == '-1\n') and (lines[i+1].lower() == 'kinetic energy\n')]
                        #blocks_4 = [i for i in range(len(lines)) if ('-1' in lines[i])]
                        print(lines[9].rstrip('\r\n'))
                        print(repr(lines[6]))
                        print(repr(lines[56]))
                        self.first_comment[(counter,subcounter)] = lines[9].rstrip('\r\n')
                        if len(blocks_4) > 1: # If trigged, i in a few lines must be redefined.
                            print('*** Interesting! More than 1 scan has been detected in above single file!')

                        # Copy data points
                        print(blocks_4)
                        i = blocks_4[0]
                        self.line_energy = float(lines[i-17]) # Will this index change if missing a transition comment?
                        print(self.line_energy)
                        self.mode[(counter,subcounter)] = lines[i-11].rstrip('\r\n')
                        self.mode_value[(counter,subcounter)] = float( lines[i-10].rstrip('\r\n') )
                        self.dwell[(counter,subcounter)] = float( lines[i+9].rstrip('\r\n') )
                        self.repeats[(counter,subcounter)] = float( lines[i+10].rstrip('\r\n') )
                        data_points = int( lines[i+16] )
                        self.cps[(counter,subcounter)] = np.zeros(data_points)
                        E_step = float( lines[i+4].rstrip('\r\n') )
                        E_start = float( lines[i+3].rstrip('\r\n') )
                        self.KE[(counter,subcounter)] = np.arange(data_points)*E_step + E_start
                        self.BE[(counter,subcounter)] = self.line_energy - self.KE[(counter,subcounter)]
                        for counter_inner in range(data_points):
                            self.cps[(counter,subcounter)][counter_inner] = float( lines[i+19+counter_inner] )/self.dwell[(counter,subcounter)]/self.repeats[(counter,subcounter)]

        # After load of data
        if self.line_energy == 0:
            print('NB: Anode material not detected! "self.BE" will contain the kinetic energy.')

def Separate_plots(data, spacer=50):
    """Take a set of data scans and return them separated along the y-axis """
    new_cps = {}

    # loop over every region
    for region in data.filekeys.keys():

        # loop over every scan in a region 
        for i in data.filekeys[region]:
            if i == 0:
                new_cps[region, i] = data.cps[region, i]
                continue
            difference = data.cps[region, i] - new_cps[region, i-1]
            minimum = min(difference)
            if minimum > spacer:
                new_cps[region, i] = data.cps[region, i]
                continue
            elif minimum > 0:
                new_cps[region, i] = data.cps[region, i] + (spacer - minimum)
                continue
            else:
                new_cps[region, i] = data.cps[region, i] + abs(minimum) + spacer
                continue
    return new_cps



