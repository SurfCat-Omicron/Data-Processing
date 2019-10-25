import numpy as np
import matplotlib.pyplot as plt
import common_toolbox as ct

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
        self.KE = XpsDict()
        self.BE = XpsDict()
        self.cps = XpsDict()
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
                print('Old VAMAS format - FIX ME!!')
                self.format = 'Old VAMAS'
                # Find identifiers
                blocks_4 = [i for i in range(len(lines)) if (lines[i].strip() == '-1') and (lines[i+1].lower().strip() == 'kinetic energy')]
                blocks_2 = [i for i in range(len(lines)) if (lines[i].strip() == 'XPS') and (lines[i+2].strip() == '1253.6')]
                print('Using old format: Mg anode hardcoded into script!')
                # Copy data points
                self.scans = len(blocks_4)
                for counter in range(self.scans):
                    i = blocks_4[counter]
                    self.mode[counter, 0] = lines[i-11]
                    self.mode_value[counter, 0] = float(lines[i-10])
                    self.dwell[counter, 0] = float(lines[i+9])
                    self.repeats[counter, 0] = float(lines[i+10])
                    data_points = int(lines[i+16])
                    self.cps[counter, 0] = np.zeros(data_points)
                    e_step = float(lines[i+4])
                    e_start = float(lines[i+3])
                    self.KE[counter, 0] = np.arange(data_points)*e_step + e_start
                    self.BE[counter, 0] = self.line_energy - self.KE[counter, 0]
                    for counter_inner in range(data_points):
                        self.cps[counter, 0][counter_inner] = float(lines[i+19+counter_inner])/self.dwell[counter, 0]/self.repeats[counter, 0]
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
                        ydata = np.zeros(data_points)
                        e_step = float(lines[i+4])
                        e_start = float(lines[i+3])
                        # Kinetic energy from file
                        self.KE[(counter, subcounter)] = np.arange(data_points)*e_step + e_start
                        # Binding energy from kinetic energy
                        self.BE[(counter, subcounter)] = self.line_energy - self.KE[(counter, subcounter)]
                        # Analyzer work function
                        Phi = 4.5
                        self.KE[(counter, subcounter)] -= Phi
                        for counter_inner in range(data_points):
                            ydata[counter_inner] = float(lines[i+19+counter_inner])/self.dwell[(counter, subcounter)]/self.repeats[(counter, subcounter)]
                        self.cps[(counter, subcounter)] = ydata

        # After load of data
        if self.line_energy == 0:
            print('NB: Anode material not detected! "self.BE" will contain the kinetic energy.')

class XpsDict(dict):
    """Custom dict to simplify indexation"""
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    def __getitem__(self, index):
        """Assure that the averaged scans are returned if an integer index
is supplied regardless of how the data was exported"""
        if isinstance(index, int):
            # The averaged data set is wanted
            data = None
            for key, value in dict.items(self):
                if key[0] == index:
                    if data is None:
                        data = value
                        counter = 1
                    else:
                        data += value
                        counter += 1
            if data is None:
                raise IndexError('Index ({}, x) doesn\'t exist.'.format(index))
            return data/float(counter)
        elif isinstance(index, tuple):
            # A single scan is specified
            return dict.__getitem__(self, index)


def plot(data, plots, plot_args, plot_kwargs, xaxis='binding energy', rc=False, y_align=None, grid=False):
    """Take dicts of options to produce standardized plots"""

    # Inputs
    if rc:
        try:
            import rcparam
        except ImportError:
            print('No module "rcparam" found. Using default plot settings.')
    if xaxis.lower().startswith('b'):
        xaxis = 'BE'
        xlabel = 'Binding energy (eV)'
        flip = True
    elif xaxis.lower().startswith('k'):
        xaxis = 'KE'
        xlabel = 'Kinetic energy (eV)'
        flip = False
    else:
        msg = 'xaxis="{}" not understood. Allowed (binding energy/kinetic energy)'
        raise ValueError(msg.format(xaxis))

    # Plots
    for plot in plots:
        plt.figure()
        plt.title(plot['title'])
        plt.xlabel(xlabel)
        plt.ylabel('CPS')
        for key, index in plot['data']:
            if xaxis == 'BE':
                x = data[key].BE[index]
            elif xaxis == 'KE':
                x = data[key].KE[index]
            y = data[key].cps[index]
            if y_align == 'right':
                edge = np.average(y[-5:])
            elif y_align == 'left':
                edge = np.average(y[:5])
            else:
                edge = 0
            y = y - edge
            plt.plot(x, y, plot_args[key])#, **plot_kwargs[key])
        if flip:
            ax = plt.gca()
            ct.flip_x(ax)
        if grid:
            print('drawing grid')
            #plt.grid(True, which='both', axis='X', color='k', linewidth=2)
            plt.grid(True, 'major', 'both')


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

class LoadSet():
    """Load a set of XPS data """

    def __init__(self, input_files=[]):
        """input_files (list of dicts/dict): 'path' to filenames"""
        if isinstance(input_files, dict):
            self.input_files = [input_files]
        elif isinstance(input_files, list):
            #if not isinstance(input_files[0], dict):
            #    raise InputError('input_files must be a list of dictionaries!')
            self.input_files = input_files

        # Data variables
        self.data = {}

    def add_data(self, input_files):
        """Add data to load."""
        if isinstance(input_files, dict):
            self.input_files.append(input_files)
        elif isinstance(input_files, list):
            if not isinstance(input_files[0], dict):
                raise ValueError('input_files must be a list of dictionaries!')
            for item in input_files:
                self.input_files.append(item)

    def load(self):
        """Compile and return dictionary of data"""

        # Load data
        loaded_files = [x.filename for x in self.data.values()]
        for input_file in self.input_files:
            path = input_file['path']
            if len(path) > 0:
                if path[-1] != '/':
                    path += '/'
            for key, filename in input_file.items():
                if key == 'path':
                    continue
                if path+filename in loaded_files:
                    continue
                else:
                    print('Input file: {}'.format(path+filename))
                self.data[key] = Experiment(path + filename)

        # Return data
        return self.data
