#import matplotlib.pyplot as plt
import numpy as np

class Experiment():
    """Load an XPS experiment exported as text or VAMAS file.

Author: Jakob Ejler Sorensen
Version: 1.2
Date: 2019 February 20
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
        f_handle = open(filename, 'r')
        print('Input file: ' + filename)
        lines = f_handle.readlines()
        f_handle.close()

        # Read data from textfile:
        if filename.endswith('.txt'):
            raise ImportError('Text file import has not been implemented yet!')
        # Read data from VAMAS block file
        elif filename.endswith('.vms'):
            # Open filename with XPS data
            f_handle = open(filename, 'r')
            lines = f_handle.readlines()
            f_handle.close()
            # Old format:
            if lines[6].lower().startswith('experiment type'):
                self.format = 'Old VAMAS'
                # Find identifiers
                blocks_4 = [i for i in range(len(lines)) if (lines[i].strip() == '-1') and (lines[i+1].lower().strip() == 'kinetic energy')]
                blocks_2 = [i for i in range(len(lines)) if (lines[i].strip() == 'XPS') and (lines[i+2].strip() == '1253.6')]
                print('Using old format: Mg anode hardcoded into script!')
                # Copy data points
                self.scans = len(blocks_4)
                for counter in range(self.scans):
                    i = blocks_4[counter]
                    self.mode[counter] = lines[i-11]
                    self.mode_value[counter] = float(lines[i-10])
                    self.dwell[counter] = float(lines[i+9])
                    self.repeats[counter] = float(lines[i+10])
                    data_points = int(lines[i+16])
                    self.cps[counter] = np.zeros(data_points)
                    e_step = float(lines[i+4])
                    e_start = float(lines[i+3])
                    self.KE[counter] = np.arange(data_points)*e_step + e_start
                    self.BE[counter] = self.line_energy - self.KE[counter]
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = float(lines[i+19+counter_inner])/self.dwell[counter]/self.repeats[counter]
            # New format
            elif lines[6].lower().startswith('created with'):
                self.format = 'New VAMAS'
                original_ending = filename.split('--')[1]
                new_ending = '--{}_{}-Detector_Region.vms'
                filen = filename.rstrip(original_ending)
                self.filename = filen

                # Detect number of regions in XPS data
                for counter in range(100):
                    for subcounter in range(1000):
                        try:
                            f_handle = open(filen + new_ending.format(counter+1, subcounter+1), 'r')
                            if not counter in self.filekeys:
                                self.filekeys[counter] = []
                            self.filekeys[counter].append(subcounter)
                            f_handle.close()
                        except IOError:
                            break
                    if subcounter == 0:
                        break

                # Open filename with XPS data
                for counter in self.filekeys.keys():
                    for subcounter in self.filekeys[counter]:
                        new_filename = filen + new_ending.format(counter+1, subcounter+1)
                        f_handle = open(new_filename, 'r')
                        lines = f_handle.readlines()
                        f_handle.close()
                        print('Loading file: ' + new_filename)
                        print('New format')
                        # ID block 4 by their "-1" line
                        blocks_4 = [i for i in range(len(lines)) if (lines[i].strip() == '-1') and (lines[i+1].lower().strip() == 'kinetic energy')]
                        print(lines[9].strip()) # Comment
                        self.first_comment[(counter, subcounter)] = lines[9].strip()
                        if len(blocks_4) > 1: # If trigged, i in a few lines must be redefined.
                            print('*** Interesting! More than 1 scan has been detected in above single file!')

                        # Copy data points
                        i = blocks_4[0]
                        self.line_energy = float(lines[i-17])
                        self.mode[(counter, subcounter)] = lines[i-11].strip()
                        self.mode_value[(counter, subcounter)] = float(lines[i-10])
                        self.dwell[(counter, subcounter)] = float(lines[i+9])
                        self.repeats[(counter, subcounter)] = float(lines[i+10])
                        data_points = int(lines[i+16])
                        self.cps[(counter, subcounter)] = np.zeros(data_points)
                        e_step = float(lines[i+4])
                        e_start = float(lines[i+3])
                        self.KE[(counter, subcounter)] = np.arange(data_points)*e_step + e_start
                        self.BE[(counter, subcounter)] = self.line_energy - self.KE[(counter, subcounter)]
                        for counter_inner in range(data_points):
                            self.cps[(counter, subcounter)][counter_inner] = float(lines[i+19+counter_inner])/self.dwell[(counter, subcounter)]/self.repeats[(counter, subcounter)]

        # After load of data
        if self.line_energy == 0:
            print('NB: Anode material not detected! "self.BE" will contain the kinetic energy.')

def separate_plots(data, spacer=50):
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
