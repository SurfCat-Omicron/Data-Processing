# pylint: disable=invalid-name,too-many-arguments
"""FIXME"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import simps
#from scipy.stats import linregress

class Experiment():
    """Load an ISS experiment exported as text or VAMAS file.

Author: Jakob Ejler Sorensen
Version: 3.0
Date: 2017 June 23
    """

    def __init__(self, filename, mass=4, theta=146.7, E0=1000):
        """Initialize the class"""
        # Constants
        self.settings = dict()
        self.settings['mass'] = mass
        self.settings['theta'] = theta
        self.settings['E0'] = E0

        # Initialize variables
        self.energy = dict()
        self.cps = dict()
        self.dwell = dict()
        self.mode = dict()
        self.mode_value = dict()
        self.filename = filename

        # Convenience function variables
        self.peak_positions = None
        self.peak_heights_raw = None
        self.peak_heights_bg = None
        self.background = None
        self.background_settings = {
            'type': None,
            'ranges': None,
            'on': False,
            }

        # Read data from textfile:
        if filename.endswith('.txt'):
            # Open filename with ISS data
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            self.format = 'Text file'

            start_points = [i for i, line in enumerate(lines) if line == 'Energy\tCounts\r\n']
            self.scans = len(start_points)
            if self.scans == 0:
                raise ImportError('File apparently empty!')

            if lines[0].lower().startswith('note'):
                self.note = lines[0].split('=')[1].lstrip(' ')
            # Copy data points
            counter = 0
            for start in start_points:
                line = lines[start-4].split('\t')
                for i, word in enumerate(line):
                    if word.lower() == 'dwell':
                        self.dwell[counter] = float(lines[start-3].split('\t')[i])
                if not start == start_points[-1]:
                    interval = range(start+1, start_points[counter+1]-4)
                else:
                    interval = range(start+1, len(lines))
                self.energy[counter] = np.zeros(len(interval))
                self.cps[counter] = np.zeros(len(interval))
                counter_inner = 0
                for index in interval:
                    line = lines[index].rstrip().split('\t')
                    self.energy[counter][counter_inner] = float(line[0])
                    self.cps[counter][counter_inner] = float(line[1])
                    counter_inner += 1
                self.cps[counter] = self.cps[counter]/self.dwell[counter]
                counter += 1
        # Read data from old VAMAS block file
        elif filename.endswith('.vms'):
            # Open filename with ISS data
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            # Old format:
            if lines[6].lower().startswith('experiment type'):
                self.format = 'Old VAMAS'
                print('Loading file: ' + filename)
                blocks_4 = [i for i, line in enumerate(lines) if (line.strip() == '-1') \
and (lines[i+1].lower().strip() == 'kinetic energy')]
                blocks_2_ISS = [i for i, line in enumerate(lines) if (line.strip() == 'ISS') \
and (lines[i+1].strip() == '')]
                print(lines[9].rstrip())
                self.scans = len(blocks_4)
                if len(blocks_4) == int(lines[9].rstrip()) \
and len(blocks_4) == len(blocks_2_ISS):
                    self.scans = len(blocks_4)
                else:
                    msg = 'Error: Identified {} "Block 4", {} "Block 2", but "Block 1" says: {}'
                    msg = msg.format(len(blocks_4), len(blocks_2_ISS), int(lines[9].rstrip()))
                    raise ImportError(msg)

                # Copy data points
                self.note = dict()
                for counter, block in enumerate(blocks_4):
                    if not len(lines[blocks_2_ISS[counter] - 1]) == 5:
                        self.note[counter] = lines[blocks_2_ISS[counter] - 1].rstrip()
                    else:
                        self.note[counter] = ''
                    self.mode[counter] = lines[block-11].rstrip()
                    self.mode_value[counter] = float(lines[block-10].rstrip())
                    self.dwell[counter] = float(lines[block+9].rstrip())
                    data_points = int(lines[block+16])
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float(lines[block+4].rstrip())
                    E_start = float(lines[block+3].rstrip())
                    self.energy[counter] = np.arange(data_points)*E_step + E_start
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = float(lines[block+19+counter_inner]) \
/self.dwell[counter]
                self.note[counter] = ''
                print(self.energy.keys())
                print('Comments: {}'.format(self.note))
                print('Dwell time: {}'.format(self.dwell))
                print('Modes: {}'.format(self.mode))
                print('Mode values: {}'.format(self.mode_value))
            # New format:
            if lines[6].lower().startswith('created with'):
                self.format = 'New VAMAS'
                ENDING = '_1-Detector_Region.vms'
                filen = filename.rstrip(ENDING)
                counter = 0
                while True:
                    try:
                        f = open(filen + '--' + str(counter+1) + ENDING)
                        counter += 1
                    except:
                        #print('{} files detected of series:'.format(counter))
                        #print('* ' + filen + '--' + str(1) + ENDING + ' *')
                        COUNTER = counter
                        break
                # Open filename with ISS data
                self.scans = COUNTER
                for counter in range(COUNTER):
                    new_filename = filen + '--' + str(counter+1) + ENDING
                    f = open(new_filename, 'r')
                    lines = f.readlines()
                    f.close()
                    print('Loading file: ' + new_filename)
                    blocks_4 = [i for i, line in enumerate(lines) if (line.rstrip() == '-1') \
and (lines[i+1].lower().rstrip() == 'kinetic energy')]
                    #print(lines[9].rstrip())
                    if len(blocks_4) > 1:
                        print('*** Interesting! More than 1 scan has been detected in above file!')
                    # Copy data points
                    i = blocks_4[0]
                    self.mode[counter] = lines[i-11].rstrip()
                    self.mode_value[counter] = float(lines[i-10].rstrip())
                    self.dwell[counter] = float(lines[i+9].rstrip())
                    data_points = int(lines[i+16])
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float(lines[i+4].rstrip())
                    E_start = float(lines[i+3].rstrip())
                    self.energy[counter] = np.arange(data_points)*E_step + E_start
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = float(lines[i+19+counter_inner]) \
/self.dwell[counter]


    def ConvertEnergy(self, mass):
        """Converts a measured energy to mass of surface atom
corresponding the settings stored in the experiment.
        """
        angle = self.settings['theta'] * np.pi/180
        return self.settings['E0'] * ((self.settings['mass']*np.cos(angle) + \
np.sqrt(mass**2 - self.settings['mass']**2*np.sin(angle)**2))/(mass + self.settings['mass']))**2


    def PlotAllScans(self, exclude=[None], color=None):
        """Plot all elements in file in single figure."""
        selection = [i for i in range(self.scans) if not i in exclude]
        if not color:
            for i in selection:
                plt.plot(self.energy[i], self.cps[i])
        else:
            for i in selection:
                plt.plot(self.energy[i], self.cps[i], color=color)
        plt.xlabel('Kinetic energy (eV)')
        plt.ylabel('Counts per second')


    def Normalize(self, interval='Total', exclude=[None], unit='Mass'):
        """Normalize to highest value in interval=[value1, value2]"""
        if isinstance(interval, int):
            self.normalization_criteria = interval
        elif isinstance(interval, str):
            if interval == 'Total':
                self.normalization_criteria = 'all'
            elif interval == 'Au':
                self.normalization_criteria = 196.
        if not isinstance(interval, list):
            if self.normalization_criteria == 'all':
                selection = [i for i in range(self.scans) if (not i in exclude)]
                for __counter in selection:
                    total = simps(self.cps[__counter], self.energy[__counter])
                    self.cps[__counter] /= total
            else:
                interval = [0, 0]
                interval[0] = self.ConvertEnergy(self.normalization_criteria) - 10
                interval[1] = self.ConvertEnergy(self.normalization_criteria) + 10
                selection = [i for i in range(self.scans) if (not i in exclude) and \
                             (not interval[0] > max(self.energy[i])) and (not interval[1] < min(self.energy[i]))]
                for __counter in selection:
                    range_1 = np.where(self.energy[__counter] < interval[1])[0]
                    range_2 = np.where(self.energy[__counter] > interval[0])[0]
                    energy_range = np.intersect1d(range_1, range_2)
                    value = max(self.cps[__counter][energy_range])
                    self.cps[__counter] = self.cps[__counter]/value


    def AddMassLines(self, masses, offset=0, color='k', labels=True):
        """Add vertical lines for mass references."""
        energies = self.ConvertEnergy(np.array(masses))
        ax = plt.gca()
        [x1, x2, y1, y2] = ax.axis()
        for energy, mass in zip(energies, masses):
            ax.axvline(x=energy-offset, ymin=0, ymax=1, linestyle='dotted', color=color)
            if labels:
                ax.text(float(energy)/x2, 0.95, 'm-{}'.format(mass), transform=ax.transAxes)


    def AddRegions(self):
        """Add regions indicating the whereabouts of 3d, 4d, 5d metals and the
lanthanides and actinides."""
        ax = plt.gca()
        d3 = [45, 65]
        d4 = [89, 112]
        d5 = [178, 201]
        lant = [139, 175]
        act = [227, 260]
        for i in [d3, d4, d5]:
            ax.axvspan(xmin=self.ConvertEnergy(i[0]), xmax=self.ConvertEnergy(i[1]),
                       color='k', alpha=0.2)
        for i in [lant, act]:
            ax.axvspan(xmin=self.ConvertEnergy(i[0]), xmax=self.ConvertEnergy(i[1]),
                       color='y', alpha=0.2)


# Convenience functions
def iss_to_time_axis(iss_data=[], ax=None):
    """Fold a series of iss measurements unto a time axis to separate them graphically"""

    colors = ['k', 'r', 'g', 'b', 'm', 'y', 'c']

    # Allow to not pass single argument as a list object
    if isinstance(iss_data, np.ndarray):
        iss_data = [iss_data]

    # Create axes object
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    # Plot iss
    separator = 0
    print('Time-plot of ISS made:\n-----')
    for j, iss in enumerate(iss_data):
        print(j, iss.filename)
        for i in iss.energy.keys():
            x, y = iss.energy[i], iss.cps[i]
            ax.plot(x + separator, y, color=colors[j])
            separator += max(x) + 10
    ax.set_xlabel('Energy (eV) -> Time')
    ax.set_ylabel('Signal (cts/s)')

    # return handles to figure
    fig = plt.gcf()
    return fig, ax

def get_range(X, limit):
    """Return index of list X where limit[0] < x < limit[1] """
    try:
        A = np.where(X > limit[0])[0][0]
        B = np.where(X < limit[1])[0][-1]
        return np.arange(A, B)
    except IndexError:
        return None

def get_peaks(iss_data=[], peaks=None):
    """FIXME"""
    AVG = 3.5
    for iss in iss_data:
        iss.peak_positions = peaks
        iss.peaks_raw = dict()
        if iss.background_settings['on']:
            iss.peaks_bg = dict()
        for key in iss.energy.keys():
            iss.peaks_raw[key] = dict()
            if iss.background_settings['on']:
                iss.peaks_bg[key] = dict()
            for element, pos in peaks.items():
                indice = get_range(iss.energy[key], [pos-AVG, pos+AVG])
                if indice is None:
                    #print('Peaks', iss.filename, key, (element, pos), 'None')
                    average = -500
                    average = np.nan
                else:
                    average = np.average(iss.cps[key][indice])
                    #print('Peaks', iss.filename, key, (element, pos), average)
                iss.peaks_raw[key][element] = average
                if iss.background_settings['on']:
                    if indice is None:
                        iss.peaks_bg[key][element] = np.nan
                    else:
                        iss.peaks_bg[key][element] = iss.peaks_raw[key][element] - \
iss.background[key][int(np.average(indice))]


class LoadSet():
    """Load a set of ISS data """
    def __init__(self, input_files=[]):
        """input_files (list of dicts/dict): 'path' to filenames"""
        if isinstance(input_files, dict):
            self.input_files = [input_files]
        elif isinstance(input_files, list):
            #if not isinstance(input_files[0], dict):
            #    raise InputError('input_files must be a list of dictionaries!')
            self.input_files = input_files

    def add_data(self, input_files):
        """Add data to load."""
        if isinstance(input_files, dict):
            self.input_files.append(input_files)
        elif isinstance(input_files, list):
            if not isinstance(input_files[0], dict):
                raise InputError('input_files must be a list of dictionaries!')
            for item in input_files:
                self.input_files.append(item)

    def load(self, normalize=None, exclude=[], unit='Mass'):
        """Compile and return dictionary of data"""
        return_dict = {}

        # Load data
        for input_file in self.input_files:
            path = input_file['path']
            if len(path) > 0:
                if path[-1] != '/':
                    path += '/'
            for key, filename in input_file.items():
                if key == 'path':
                    continue
                return_dict[key] = Experiment(path + filename)

        # Normalize data if chosen
        if not normalize is None:
            for key in return_dict:
                return_dict[key].Normalize(normalize, exclude, unit)

        # Return data
        return return_dict

def plot(data, index, ax=None, args=[], kwargs={}):
    """Plot index from Experiment 'data'"""
    if ax is None:
        ax = plt.gca()
    x, y = data.energy[index], data.cps[index]
    ax.plot(x, y, *args, **kwargs)

def subtract_backgrounds(iss_data=[], ranges=[], btype='linear', avg=3):
    """Subtract a linear background from defined 'ranges'.
Return data above backgrounds."""

    AVG = avg + 0.5
    for iss in iss_data:
        iss.background = dict()
        for key in iss.energy.keys():
            background = np.zeros(len(iss.cps[key]))
            background[:] = iss.cps[key]
            for limit in ranges:
                indice = get_range(iss.energy[key], limit)
                # if first index is chosen
                # OR
                # if last ten indice are included
                #if indice[0] < limit[0] or indice[-1] > limit[1]:
                if indice[0] == 0 or indice[-1] > len(iss.energy[key]) - 10:
                    #print(len(iss.energy[key]))
                    print('Uhh', indice[0], indice[-1], limit)
                    background[indice] = 0
                    #print('BG **', iss.filename, key, indice)
                elif len(indice) == 0:
                    print('BG *', iss.filename, key, indice)
                else:
                    y1 = np.average(iss.cps[key][indice[0]-AVG:indice[0]+AVG])
                    y2 = np.average(iss.cps[key][indice[-1]-AVG:indice[-1]+AVG])
                    a_coeff = (y2-y1)/(limit[1]-limit[0])
                    b_coeff = y1 - a_coeff*limit[0]
                    print(y1, y2, a_coeff, b_coeff)
                    background[indice] = iss.energy[key][indice]*a_coeff + b_coeff
            iss.background[key] = background

        # Tell settings that function is done
        iss.background_settings['type'] = btype
        iss.background_settings['ranges'] = ranges
        iss.background_settings['on'] = True

def reorganize_peaks(iss_data=[], bg=True):
    """Take the list of dicts() of peaks and gather into a dict of elements
TODO: add possibility of only selecting specific sets (iss.energy[i])"""
    peaks = dict()
    for element in iss_data[0].peak_positions.keys():
        if bg:
            peaks[element] = [iss.peaks_bg[i][element] for iss in iss_data \
for i in iss.peaks_bg.keys()]
            peaks[element] = np.array(peaks[element])
        else:
            peaks[element] = [iss.peaks_raw[i][element] for iss in iss_data \
for i in iss.peaks_raw.keys()]
            peaks[element] = np.array(peaks[element])
    total = np.zeros(len(peaks[element]))
    for element in peaks.keys():
        total += peaks[element]
    return peaks, total

def convolute_peaks(iss_data=[], peaks=[], sigma=[30, 30], amp=[1000, 1000]):
    """
Take a set of iss data with background subtraction and fit a dictionary of peak positions.
Return: None. Peaks are saved as attributes of ISS object under 'peak_fit' as a dict.
"""
    def gauss(x, x0, sigma, amp):
        """Gaussian function"""
        return amp * np.exp(-((x - x0)**2)/(2*sigma**2))

    def fit_gauss(x, y, x01, x02, sigma1, sigma2, amp1, amp2):

        def double_gauss(x, sigma1, sigma2, amp1, amp2):
            return gauss(x, x01, sigma1, amp1) + gauss(x, x02, sigma2, amp2)

        param = [sigma1, sigma2, amp1, amp2]
        ret, _null = curve_fit(double_gauss, x, y, p0=param)
        return ret

    for iss in iss_data:
        iss.peak_fit = dict()
        for i in range(iss.scans):
            signal = iss.cps[i] - iss.background[i]
            energy = iss.energy[i][:]
            iss.peak_fit[i] = fit_gauss(energy, signal, peaks[0], peaks[1], sigma[0], sigma[1], \
amp[0], amp[1])
            print(iss.peak_fit[i])
