# pylint: disable=
import numpy as np

def smooth(data, width=1):
    """Average `data` with `width` neighbors"""
    smoothed_data = np.zeros(len(data))
    smoothed_data[width:-width] = data[2*width:]
    for i in range(2*width):
        smoothed_data[width:-width] += data[i:-2*width+i]
        if i < width:
            smoothed_data[i] = sum(data[0:i+width+1])/len(data[0:i+width+1])
            smoothed_data[-1-i] = sum(data[-1-i-width:])/len(data[-1-i-width:])
    smoothed_data[width:-1-width] = smoothed_data[width:-1-width]/(2*width+1)
    return smoothed_data

def get_averaged_gradient(current, width=2):
    """Average gradient over `width` neighbors."""
    """Interpretation:
    gradient[i] related to time[i] is the change at time[i] to time[i+1]
    gradient[i] related to time[i+1] is the change at time[i+1] from time[i]
    averaged_gradient is the gradient averaged over `width` neighbors
    """

    gradient = current[1:] - current[0:-1]
    averaged_gradient = np.zeros(len(current))
    # calculate middle
    denominator = 2*width + 1
    for i in range(denominator-1):
        averaged_gradient[width:-width-1] += gradient[i:i-2*width]/denominator
    averaged_gradient[width:-width-1] += gradient[i+1:]/denominator
    # append ends
    averaged_gradient[0:width] = averaged_gradient[width]
    averaged_gradient[-width-1:] = averaged_gradient[-width-2]
    return averaged_gradient

def parse_bool_string(string): # Interpret string as true or false
    """Convert string to boolean value."""
    return (string[0].upper() == 'T') or (string.upper() == 'ON') or (string[0].upper() == 'Y')

def to_pmol(number):
    """Convert a number to [pmol] unit."""
    return number*1e12/6.022e23

class IntegrateCurrent(object): # pylint: disable=useless-object-inheritance
    """Integrate current in an 'omicron' deposition measurement to get total charge and coverage.

    Parameters
    ----------
    data : ndarray
        time as first column in `data`
        current as second column in `data`
    parameter_string : string
        semicolon-separated 'key=value" pairs:

        MODEL -- model to calculate coverage from (NP/SA)
        SA_DENSITY -- "SA" single atom density (unit: cm**-2)
        PARTICLE_DIAMETER -- "NP" diameter of nanoparticles (unit: nm)
        APERTURE_DIAMETER -- "NP/SA" diameter of deposition aperture/area (unit: mm)
        TARGET -- target coverage used for time estimate (unit: percent, default: 0)
        FIRST_LIMIT -- approximate amplitude of first registered current (unit: pA)
        SENSITIVITY_LIMIT -- fraction of deposition current used to detect new background current
        SENSITIVITY_FILTER -- factor for noise filter. A high factor means keeping more data.
            A factor of '1' is correlated directly to the standard deviation in the beginning of
            the measurement.
        TIME -- time interval to integrate. Leave empty for full interval (unit: s)
        PLOT -- determine whether a figure should be produced (True/False)
        DEBUG -- plot a figure to help determine `SENSITIVITY_X` (True/False)
    """

    def __init__(self, parameter_string, data, plot=True):
        """Assign parameters as attributes."""
        string = parameter_string.strip(';').split(';')
        self.debugging = False
        self.model, self.target, self.time = None, None, list()
        if plot:
            import matplotlib.pyplot as plt
            self.plt = plt
            self.plot = plot
        for item in string:
            param, value = item.split('=')
            if param == 'SENSITIVITY_FILTER':
                self.sensitivity_filter = float(value)
            elif param == 'SENSITIVITY_LIMIT':
                self.sensitivity_limit = float(value)
            elif param == 'FIRST_LIMIT':
                self.first_limit = float(value)*1e-12
            elif param == 'PARTICLE_DIAMETER':
                radius_particle = float(value)/2.*1e-9
                self.area_particle = np.pi * ((radius_particle)**2)
            elif param == 'APERTURE_DIAMETER':
                radius_aperture = float(value)/2.*1e-3
                self.area_aperture = np.pi * ((radius_aperture)**2)
            elif param == 'DEBUG':
                self.debugging = parse_bool_string(value)
            elif param == 'MODEL':
                if value in ['NP', 'SA']:
                    self.model = value
                else:
                    print('MODEL must have either "NP" or "SA"')
            elif param == 'SA_DENSITY':
                self.sa_density = float(value)
            elif param == 'TARGET':
                self.target = [float(value)]
                if not self.target[0] > 0:
                    self.target = np.array([5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300])
            elif param == 'TIME':
                value = value.strip('[]')
                if value:
                    self.time = [float(val) for val in value.split(',')]
            else:
                print('par/val pair ({0}/{1}) not recognized'.format(param, value))

        # Store data
        self.data = data
        self.dep_rate = None
        self.present_coverage = None

    def convert_coverage_to_number(self, coverage):
        """Convert target coverage from percent to number of particles."""
        if self.model is None:
            print('Model not chosen!')
            number = None
        elif self.model == 'SA':
            print('Single atoms model not implemented yet!')
            number = None
        elif self.model == 'NP':
            number = coverage/100.*self.area_aperture/self.area_particle
        return number

    def convert_number_to_coverage(self, number):
        """Convert coverage from number of particles to percent."""
        if self.model is None:
            print('Model not chosen!')
            coverage = None
        elif self.model == 'SA':
            print('Single atoms model not implemented yet!')
            coverage = None
        elif self.model == 'NP':
            coverage = number*100/self.area_aperture*self.area_particle
        return coverage

    def integrate(self):
        """Integrate the current and return the number of particles."""

        # Filter data
        self.filter_data()

        # Separate data in regions
        leak_current = self.separate_data()

        # Integrate deposition current
        e = 1.602e-19 # pylint: disable=invalid-name
        net_current = np.abs(self.data[:, 1] - leak_current[:, 1])
        integral = np.sum(net_current[1:] * np.diff(self.data[:, 0]))/e
        dep_rate = net_current[np.where(net_current > 1e-13)][-100:]
        dep_rate = np.average(dep_rate)/e

        # Estimate time to next coverage step
        present_coverage = self.convert_number_to_coverage(integral)
        for coverage in self.target:
            if coverage >= present_coverage:
                target_coverage = coverage
                break
            else:
                target_coverage = self.target[-1]
        target_number = self.convert_coverage_to_number(target_coverage)

        # If specific target chosen or over 300 %, else...
        if target_number < integral:
            msg = 'Target coverage ({0} %) already exceeded: {1} %.'
            print(msg.format(target_coverage, present_coverage))
            print('Number of clusters: {0} pmol.'.format(to_pmol(integral)))
        else:
            msg = 'Total charge deposited: {} pmol.\n'.format(to_pmol(integral))
            msg += 'Present coverage: {} %\n'.format(present_coverage)
            msg += 'Time until {} % coverage: {} hr {} min {} sec'
            remaining = target_number - integral
            time_left = remaining/dep_rate
            time_hr = int(time_left/3600)
            time_min = int((time_left - time_hr*3600)/60)
            time_sec = int(time_left - time_hr*3600 - time_min*60)
            print(msg.format(target_coverage, time_hr, time_min, time_sec))

    def separate_data(self):
        """Separate data into regions of background and deposition current."""

        # Identify type of first region
        counter = 0
        offset = 1
        limit = self.first_limit*self.sensitivity_limit
        try:
            while True:
                index = offset + counter
                if self.data[index+1, 1] - self.data[index, 1] > limit:
                    first_region = 'deposition'
                    break
                elif self.data[index+1, 1] - self.data[index, 1] < -limit:
                    first_region = 'background'
                    break
                counter += 1
            msg = 'First region detected: "{0}" at filtered index {1} = {2} seconds'
            print(msg.format(first_region, index, round(self.data[index, 0], 2)))
        except IndexError:
            print(' *** Index exceeded when searching for first region. ***')

        # Find remaining regions
        deposition, leak = dict(), dict()
        counter_current, counter_leak = 0, 0
        if first_region == 'background':
            depo = False
        else:
            depo = True
        previous = -1
        for i, gradient in enumerate(self.data[:-1, 1] - self.data[1:, 1]):
            # Detect end of leak measurement
            #if counter_current + counter_leak > 0:
            #    pass
            #elif limit < gradient and not depo:
            if limit < gradient and not depo:
                leak[counter_leak] = np.arange(previous+1, i+1)
                previous = i
                counter_leak += 1
                depo = True
            elif -limit > gradient and depo:
                deposition[counter_current] = np.arange(previous+1, i+1)
                previous = i
                counter_current += 1
                depo = False
                limit = self.renew_limit(i)
        # End of measurement
        if depo:
            last_region = 'deposition'
            deposition[counter_current] = np.arange(previous+1, len(self.data))
        else:
            last_region = 'background'
            leak[counter_leak] = np.arange(previous+1, len(self.data))
        msg = 'Regions detected:\n\t{0} deposition currents and\n\t{1} leak currents'
        print(msg.format(len(deposition), len(leak)))

        # Extrapolate between regions
        counter_current = len(deposition)
        extrapolated_leak = dict()
        if first_region == 'deposition':
            for i in np.arange(counter_current):
                if i == 0:
                    slope, intercept = 0, self.data[leak[i], 1].mean()
                elif (i == counter_current - 1) and (last_region == 'deposition'):
                    slope, intercept = 0, self.data[leak[i-1], 1].mean()
                else:
                    slope = (self.data[leak[i], 1].mean() - self.data[leak[i-1], 1].mean())
                    slope = slope/(self.data[leak[i][-1], 0] - self.data[leak[i-1][0], 0])
                    intercept = (self.data[leak[i], 1].mean() - slope*self.data[leak[i][-1], 0])
                extrapolated_leak[i] = slope*self.data[deposition[i], 0] + intercept
        else:
            for i in np.arange(counter_current):
                if (i == counter_current - 1) and (last_region == 'deposition'):
                    slope, intercept = 0, self.data[leak[i], 1].mean()
                else:
                    slope = (self.data[leak[i], 1].mean() - self.data[leak[i+1], 1].mean())
                    slope = slope/(self.data[leak[i][-1], 0] - self.data[leak[i+1][0], 0])
                    intercept = (self.data[leak[i], 1].mean() - slope*self.data[leak[i][-1], 0])
                extrapolated_leak[i] = slope*self.data[deposition[i], 0] + intercept
        leak_current = self.data.copy()
        for key in deposition:
            leak_current[deposition[key], 1] = extrapolated_leak[key]

        # Plot the different regions
        if self.plot:
            for key in leak.keys():
                self.plt.plot(self.data[:, 0][leak[key]], self.data[:, 1][leak[key]]/1e-12, 'go', markersize=4)
            for key in deposition.keys():
                self.plt.plot(self.data[:, 0][deposition[key]], self.data[:, 1][deposition[key]]/1e-12, 'bo', markersize=4)
            for key in deposition.keys():
                mask = deposition[key]
                self.plt.fill_between(self.data[:, 0][mask], leak_current[:, 1][mask]/1e-12, self.data[:, 1][mask]/1e-12, color='b', alpha=0.3)
        return leak_current

    def renew_limit(self, index):
        """Update the magnitude of deposition current for detecting next background."""
        left = np.average(self.data[index-5:index+1, 1])
        right = np.average(self.data[index+1:index+7, 1])
        return abs(left-right)*self.sensitivity_limit

    def filter_data(self):
        """Apply various filters to raw data."""
        # Remove any overflow data (1E+38)
        self.data = self.data[np.where(self.data[:, 1] < 1)]

        # Plot raw data
        if self.plot:
            self.plt.plot(self.data[:, 0], self.data[:, 1]/1e-12, 'ro-', markersize=2) ###
            self.plt.xlabel('Time (s)')
            self.plt.ylabel('Current (pA)')
        
        # Select a portion of the data
        if self.time:
            self.data = self.data[np.where(self.data[:, 0] >= self.time[0])]
            self.data = self.data[np.where(self.data[:, 0] <= self.time[1])]

        # Get fluctuation information
        averaged_gradient = get_averaged_gradient(self.data[:, 1], width=2)
        std = self.get_std(averaged_gradient)

        # Apply noise filter
        self.data = self.data[np.where(abs(averaged_gradient) < std*self.sensitivity_filter)]

    def get_std(self, averaged_gradient):
        """Not a disease, but the standard deviation of the measurement fluctuations."""
        counter = 0
        offset = 15
        while True:
            gradient = abs(self.data[1:offset+counter, 1]-self.data[0:offset+counter-1, 1])
            std_gradient = gradient.std()
            counter += 1
            if abs(averaged_gradient[offset+counter]) > std_gradient:
                print('Gradient computed up to index {}'.format(offset+counter))
                break
        return std_gradient

################################################################
### MAIN ###
################################################################
if __name__ == '__main__':

    #import rcparam
    from cinfdata import Cinfdata
    db = Cinfdata('omicron', use_caching=False) # pylint: disable=invalid-name
    ID = 15372
    STRING = '\
TARGET=5.0;\
MODEL=NP;\
SA_DENSITY=1;\
PARTICLE_DIAMETER=5.0;\
APERTURE_DIAMETER=4.5;\
FIRST_LIMIT=10.8;\
SENSITIVITY_LIMIT=1.;SENSITIVITY_FILTER=1.;\
TIME=[]'
    try:
        SESSION = IntegrateCurrent(STRING, db.get_data(ID))
        SESSION.integrate()
        if SESSION.plot:
            SESSION.plt.show()
    except:
        print('***\nSomething is wrong: Check input parameters or try debugging mode!!\n***')
        SESSION.plt.show()
        raise
    # for 9x9 raster pattern: ap_dia ~ 12.4 mm (120.8 mm2 ~ 11x11 mm)
    # for 5x5 raster pattern: ap_dia ~ 6.7 mm (35.3 mm2)
    #                 [*** Based on simul. 12/12-18 use ap_dia ~ 9.0mm]
    # for localized_Z pattern: ap_dia ~ 4.81 mm (18.2 mm2 ~ 5.2x3.5 mm)
