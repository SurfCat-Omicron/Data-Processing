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
        self.line_energy = 1486.6

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
                #ENDING = '_1-Detector_Region.vms'
                #filen = filename.rstrip(ENDING)
                ## Detect number of regions in XPS data
                #counter = 0
                #while True:
                    #try:
                    #    f = open(filen + '--' + str(counter+1) + ENDING)
                    #    counter += 1
                #    except:
                #        #print('{} files detected of series:'.format(counter))
                #        #print('* ' + filen + '--' + str(1) + ENDING + ' *')
                #        COUNTER = counter
                #        break
                # Open filename with XPS data
                #self.scans = COUNTER
                COUNTER = 1
                for counter in range(COUNTER):
                    #new_filename = filen + '--' + str(counter+1) + ENDING
                    new_filename = filename
                    f = open(new_filename, 'r')
                    lines = f.readlines()
                    f.close()
                    print('Loading file: ' + new_filename)
                    # ID block 4 by their "-1" line
                    blocks_4 = [i for i in range(len(lines)) if (lines[i] == '-1\n') and (lines[i+1].lower() == 'kinetic energy\n')]
                    print(lines[9].rstrip('\n'))
                    if len(blocks_4) > 1: # If trigged, i in a few lines must be redefined.
                        print('*** Interesting! More than 1 scan has been detected in above single file!')

                    # Copy data points
                    i = blocks_4[0]
                    self.mode[counter] = lines[i-11].rstrip('\r\n')
                    self.mode_value[counter] = float( lines[i-10].rstrip('\r\n') )
                    self.dwell[counter] = float( lines[i+9].rstrip('\r\n') )
                    self.repeats[counter] = float( lines[i+10].rstrip('\r\n') )
                    data_points = int( lines[i+16] )
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float( lines[i+4].rstrip('\r\n') )
                    E_start = float( lines[i+3].rstrip('\r\n') )
                    self.KE[counter] = np.arange(data_points)*E_step + E_start
                    self.BE[counter] = self.line_energy - self.KE[counter]
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = float( lines[i+19+counter_inner] )/self.dwell[counter]/self.repeats[counter]
