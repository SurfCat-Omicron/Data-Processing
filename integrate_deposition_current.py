# pylint: disable=line-too-long,missing-docstring
import sys
import matplotlib.pyplot as plt
import numpy as np
from cinfdata import Cinfdata
#########################
### FORMAT FOR PLUGIN ###
################################################################
#
#   script must be a function with two inputs:
#
#       def my_plugin(ID, string):
#
#   ID :        plot/data id number from surfcatdata
#   string :    string containing any input variables such as
#               target diameter of NPs and aperture size
#
###########################################################
### FUNCTIONS ###
###########################################################
TARGET_COVERAGE = None
if len(sys.argv) > 1:
    TARGET_COVERAGE = float(sys.argv[1])
print(TARGET_COVERAGE)
def run_script(identifier, string):
    # Get parameters --------------------------------------
    sensitivity_filter, sensitivity_limit, first_limit, plotting, debugging, radius_particle, radius_aperture = get_parameters(string)
    # Load raw data and apply filters ---------------------
    time, current = get_filtered_data(identifier, sensitivity_filter, plotting, debugging)
    # Analyze start ---------------------------------------
    limit = first_limit
    first_region = get_first_region(time, current, limit, sensitivity_limit, plotting, debugging)
    # Extrapolate leak current measurements ---------------
    extrapolated_leak, deposition, leak = extrapolate_leak(time, current, first_region, limit, sensitivity_limit, debugging)
    # Calaculate coverage and estimated remaining time ----
    coverage = get_info(time, current, deposition, extrapolated_leak, radius_particle, radius_aperture)
    # Plots -----------------------------------------------
    if plotting:
        plt.figure(1)
        for i in np.arange(len(deposition.keys())):
            plt.plot(time[deposition[i]]/60, current[deposition[i]]*1e12, 'bo')
            plt.plot(time[deposition[i]]/60, extrapolated_leak[i]*1e12, 'yo')
        for i in np.arange(leak.keys().__len__()):
            plt.plot(time[leak[i]]/60., current[leak[i]]*1e12, 'go')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Deposition current (pA)')
        plt.title('Coverage: {0:.3f}'.format(coverage))
        plt.show()
    elif debugging:
        plt.show()
###########################################################
def get_parameters(string):
    string = string.split(';') # separate string into "parameter=value" pairs
    plotting, debugging = False, False
    for i in string:
        param, value = i.split('=')
        if param == 'SENSITIVITY_FILTER':
            sensitivity_filter = float(value)
        elif param == 'SENSITIVITY_LIMIT':
            sensitivity_limit = float(value)
        elif param == 'FIRST_LIMIT':
            first_limit = float(value)*1e-12
        elif param == 'PLOT':
            plotting = parse_bool_string(value)
        elif param == 'PARTICLE_DIAMETER':
            radius_particle = float(value)/2.*1e-9
        elif param == 'APERTURE_DIAMETER':
            radius_aperture = float(value)/2.*1e-3
        elif param == 'DEBUG':
            debugging = parse_bool_string(value)
    return sensitivity_filter, sensitivity_limit, first_limit, plotting, debugging, radius_particle, radius_aperture
#----------------------------------------------------------
def parse_bool_string(string): # Interpret string as true or false
    return (string[0].upper() == 'T') or (string[0].upper() == 'ON') or (string[0].upper() == 'Y')
###########################################################
def get_filtered_data(identifier, sensitivity_filter, plotting, debugging):
    # Load raw data
    rtime, rcurrent = get_raw_data(identifier)
    #custom_filter = np.where(rtime >= 600)
    #rtime, rcurrent = rtime[custom_filter], rcurrent[custom_filter]
    if plotting:
        plt.figure(1)
        plt.plot(rtime/60, rcurrent*1e12, 'ro-')
    # Compute averaged gradient
    averaged_gradient = get_averaged_gradient(rcurrent, width=2)
    std = get_gradient(rtime, rcurrent, debugging)
    # Apply noise filter
    noise_filter = filter_noise(smooth(abs(averaged_gradient), 1), std, sensitivity_filter)
    time = rtime[noise_filter]
    current = rcurrent[noise_filter]
    if plotting:
        plt.figure(1)
        plt.plot(time/60, current*1e12, 'bo')
    if debugging:
        plt.figure(3)
        plt.title('DEBUGGING: SENSITIVITY_FILTER')
        plt.plot(rtime/60, abs(averaged_gradient)/std, 'ro-')
        plt.plot(rtime/60, np.ones(len(rtime)), 'b-')
        plt.plot(rtime/60, np.ones(len(rtime))*sensitivity_filter, 'b', linestyle='dashed')
    return time, current
#----------------------------------------------------------
def get_raw_data(identifier):
    db = Cinfdata('omicron', use_caching=False) # pylint: disable=invalid-name
    middata = db.get_data(identifier)
    metadata = db.get_metadata(identifier)
    print('\nLoading data from: {}'.format(metadata['Comment']))
    time = middata[:, 0]
    current = middata[:, 1]
    overflow_filter = filter_overflow(current)
    return time[overflow_filter], current[overflow_filter]
#----------------------------------------------------------
def filter_overflow(current): # Remove overflowed data where I > 0
    overflow_filter = np.where(current < 1)
    return overflow_filter
#----------------------------------------------------------
def get_averaged_gradient(current, width=2):
    ### INTERPRETATION: ###
    # gradient[i] related to time[i] is the change at time[i] to time[i+1]
    # gradient[i] related to time[i+1] is the change at time[i+1] from time[i]
    # averaged_gradient is the gradient averaged over *width* neighbors
    gradient = current[1:] - current[0:-1]
    num = len(current)
    averaged_gradient = np.zeros(num)
    # calculate middle
    denominator = 2*width + 1
    for i in range(denominator-1):
        averaged_gradient[width:-width-1] += gradient[i:i-2*width]/denominator
    averaged_gradient[width:-width-1] += gradient[i+1:]/denominator
    # append ends
    averaged_gradient[0:width] = averaged_gradient[width]
    averaged_gradient[-width-1:] = averaged_gradient[-width-2]
    return averaged_gradient
#----------------------------------------------------------
def get_gradient(time, current, debugging):
    counter = 0
    offset = 15
    while True:
        gradient = abs(current[1:offset+counter]-current[0:offset+counter-1])
        std_gradient = gradient.std()
        averaged_gradient = get_averaged_gradient(current[0:offset+counter+6])
        counter += 1
        if abs(averaged_gradient[offset+counter]) > std_gradient:
            print('Gradient computed up to index {}'.format(offset+counter))
            break
    if debugging:
        counter = 0
        while True:
            gradient = current[1:offset+counter]-current[0:offset+counter-1]
            new_std_gradient = gradient.std()
            averaged_gradient = get_averaged_gradient(current[0:offset+counter+6])
            plt.figure(1)
            plt.plot(time[counter+offset]/60, new_std_gradient/std_gradient, 'bo')
            counter += 1
            if abs(averaged_gradient[offset+counter]) > std_gradient:
                break
    return std_gradient
#----------------------------------------------------------
def smooth(data, width): # Smoothing function (averaging with neighbors)
    smoothed_data = np.zeros(len(data))
    smoothed_data[width:-width] = data[2*width:]
    for i in range(2*width):
        smoothed_data[width:-width] += data[i:-2*width+i]
        if i < width:
            smoothed_data[i] = sum(data[0:i+width+1])/len(data[0:i+width+1])
            smoothed_data[-1-i] = sum(data[-1-i-width:])/len(data[-1-i-width:])
    smoothed_data[width:-1-width] = smoothed_data[width:-1-width]/(2*width+1)
    return smoothed_data
#----------------------------------------------------------
def filter_noise(reduced_gradient, std_gradient, sensitivity_filter): # Filter noise
    return np.where(abs(reduced_gradient) < std_gradient*sensitivity_filter)
###########################################################
def get_first_region(time, data, first_limit, sensitivity_limit, plotting, debugging):
    counter = 0
    offset = 1
    try:
        while True:
            index = offset + counter
            if data[index+1] - data[index] > first_limit*sensitivity_limit:
                first_region = 'deposition'
                break
            elif data[index+1] - data[index] < -first_limit*sensitivity_limit:
                first_region = 'background'
                break
            counter += 1
        print('First region detected: "{0}" at filtered index {1} = {2} minutes'.format(first_region, index, round(time[index]/60, 2)))
        plot_first_limit(time, data, 0, index, first_region, first_limit, sensitivity_limit, plotting, debugging)
    except IndexError:
        print('   *** INDEXERROR: try changing FIRST_LIMIT or SENSITIVITY_LIMIT   *** ')
        plot_first_limit(time, data, 0, -2, 'background', first_limit, sensitivity_limit, plotting, debugging)
        plot_first_limit(time, data, 0, -2, 'deposition', first_limit, sensitivity_limit, plotting, debugging)
        plt.show()
    return first_region
#----------------------------------------------------------
def plot_first_limit(time, data, index, index_x, first_region, first_limit, sensitivity_limit, plotting, debugging):
    if plotting:
        plt.figure(1)
        limit = first_limit*1e12
        x = np.array([0, time[index_x+1], time[index_x+1], 0])/60
        if first_region == 'deposition':
            y = np.array([0, 0, limit, limit])+data[index]*1e12
            plt.fill(x, y, color='k', alpha=0.2)
            limit = limit*sensitivity_limit
            y = np.array([0, 0, limit, limit])+data[index]*1e12
            plt.fill(x, y, color='k', alpha=0.5)
        elif first_region == 'background':
            y = np.array([0, 0, -limit, -limit])+data[index]*1e12
            plt.fill(x, y, color='k', alpha=0.2)
            limit = limit*sensitivity_limit
            y = np.array([0, 0, -limit, -limit])+data[index]*1e12
            plt.fill(x, y, color='k', alpha=0.5)
###########################################################
def extrapolate_leak(time, data, first_region, limit, sensitivity_limit, debugging):
    deposition, leak, last_region = get_regions(time, data, first_region, limit, sensitivity_limit, debugging)
    counter_current = deposition.keys().__len__()
    extrapolated_leak = dict()
    if first_region == 'deposition':
        for i in np.arange(counter_current):
            if i == 0:
                slope, intercept = 0, data[leak[i]].mean()
            elif (i == counter_current - 1) and (last_region == 'deposition'):
                slope, intercept = 0, data[leak[i-1]].mean()
            else:
                slope = (data[leak[i]].mean() - data[leak[i-1]].mean())/(time[leak[i][-1]]-time[leak[i-1][0]])
                intercept = (data[leak[i]].mean() - slope*time[leak[i][-1]])
            extrapolated_leak[i] = slope*time[deposition[i]] + intercept
    else:
        for i in np.arange(counter_current):
            if (i == counter_current - 1) and (last_region == 'deposition'):
                slope, intercept = 0, data[leak[i]].mean()
            else:
                slope = (data[leak[i]].mean() - data[leak[i+1]].mean())/(time[leak[i][-1]]-time[leak[i+1][0]])
                intercept = (data[leak[i]].mean() - slope*time[leak[i][-1]])
            extrapolated_leak[i] = slope*time[deposition[i]] + intercept
    return extrapolated_leak, deposition, leak
#----------------------------------------------------------
def get_regions(time, data, first_region, limit, sensitivity_limit, debugging):
    # Automatically separate deposition current into regions
    # of background and actual deposition current
    deposition, leak = dict(), dict()
    counter_current, counter_leak = 0, 0
    if first_region == 'background':
        depo = False
    else:
        depo = True
    previous = -1
    if debugging:
        plt.figure(2)
        plt.title('DEBUGGING: SENSITIVITY_LIMIT')
        plt.plot(time[1:]/60, abs(data[1:]-data[0:-1]), 'bo-')
        plt.plot(time/60, np.ones(len(time))*limit, 'r-')
    # Index i represents last index of region
    for i in np.arange(len(data)-1):
        # Detect end of leak measurement
        if ((data[i+1] < data[i] - sensitivity_limit*limit) and not depo):
            leak[counter_leak] = np.arange(previous+1, i+1)
            previous = i
            counter_leak += 1
            depo = True
            if debugging and (counter_leak == 1):
                plt.plot(time[i:]/60, np.ones(len(time[i:]))*limit*sensitivity_limit, 'r', linestyle='dashed')
        # Detect end of current measurement
        elif ((data[i+1] > data[i] + sensitivity_limit*limit) and depo):
            deposition[counter_current] = np.arange(previous+1, i+1)
            previous = i
            counter_current += 1
            depo = False
            limit = renew_limit(data, i)
            if debugging:
                plt.plot(time[i:]/60, np.ones(len(time[i:]))*sensitivity_limit*limit, 'r', linestyle='dashed')
    # End of measurement
    if depo:
        last_region = 'deposition'
        deposition[counter_current] = np.arange(previous+1, len(data))
    else:
        last_region = 'background'
        leak[counter_leak] = np.arange(previous+1, len(data))
    print('\nRegions detected:\n\t{} deposition currents and\n\t{} leak currents\n'.format(len(deposition.keys()), len(leak.keys())))
    return deposition, leak, last_region
#----------------------------------------------------------
def renew_limit(data, index):
    left = data[index-5:index+1]
    left = sum(left)/len(left)
    right = data[index+1:index+7]
    right = sum(right)/len(right)
    return abs(left-right)
###########################################################
def get_info(time, current, deposition, extrapolated_leak, radius_particle, radius_aperture):
    number_of_charges, charges_per_second = integrate_current(time, current, deposition, extrapolated_leak)
    target, coverage, time_remaining_minutes, time_remaining_seconds = calculate_coverage(number_of_charges, charges_per_second, radius_particle, radius_aperture)
    #print('Deposition time: {} hr'.format(number_of_charges/charges_per_second/3600))
    print('Calculated coverage:\n   ***   {0} percent   ***\n'.format(round(coverage, 3)))
    print('Estimated remaining time to {0} percent coverage:\n   ***   {1} minutes and {2} seconds   ***'.format(target, time_remaining_minutes, round(time_remaining_seconds, 1)))
    return coverage
#----------------------------------------------------------
def integrate_current(time, current, deposition, extrapolated_leak):
    # Return number of nanoparticles
    num = len(deposition.keys())
    integration = 0
    total_time = 0
    for i in np.arange(num):
        actual_current = extrapolated_leak[i] - current[deposition[i]]
        actual_current = (actual_current[0:-1] + actual_current[1:])/2.
        integration += sum(actual_current*np.diff(time[deposition[i]]))
        total_time += sum(np.diff(time[deposition[i]]))
    e = 1.602e-19
    print('Number of charges presently measured: {0:.4} pmol'.format(integration/e/6.022e23*1e12))
    print('Time: {} hrs'.format(total_time/3600))
    #print('Charges per second: {}'.format(actual_current.mean()/e))
    return integration/e, actual_current.mean()/e
#----------------------------------------------------------
def calculate_coverage(number_of_charges, charges_per_second, radius_particle, radius_aperture):
    area_particle = np.pi * ((radius_particle)**2)
    projected_area = area_particle * number_of_charges
    area_aperture = np.pi * ((radius_aperture)**2)
    coverage = projected_area/area_aperture * 100
    if not TARGET_COVERAGE is None:
        targets = [TARGET_COVERAGE]
    else:
        targets = []
    for target in targets + [5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300]:
        if coverage < target:
            total_number_of_charges = target/100.*area_aperture/area_particle
            #relative = charges_per_second*0.128
            #print('Single step percent of current measurement: {0:.4}'.format(relative/number_of_charges*100))
            #print('Single step percent of target measurement: {0:.4}'.format(relative/total_number_of_charges*100))
            time_remaining = (total_number_of_charges - number_of_charges)/charges_per_second
            time_remaining_minutes = int(time_remaining/60)
            time_remaining_seconds = time_remaining - time_remaining_minutes*60
            break
    return target, coverage, time_remaining_minutes, time_remaining_seconds
################################################################
### MAIN ###
################################################################
if __name__ == '__main__':
    try:
        run_script(14401, 'PARTICLE_DIAMETER=5;APERTURE_DIAMETER=4.5;FIRST_LIMIT=14.0;SENSITIVITY_LIMIT=0.7;SENSITIVITY_FILTER=10.;PLOT=True;DEBUG=rue')
    except:
        print('***\nSomething is wrong: Check input parameters or try debugging mode!!\n***')
        plt.show()
        raise
    print('--- END OF SCRIPT ---')
    # for 9x9 raster pattern: ap_dia ~ 12.4 mm (120.8 mm2 ~ 11x11 mm)
    # for 5x5 raster pattern: ap_dia ~ 6.7 mm (35.3 mm2) [*** Based on simul. 12/12-18 use ap_dia ~ 9.0mm]
    # for localized_Z pattern: ap_dia ~ 4.81 mm (18.2 mm2 ~ 5.2x3.5 mm)
# DH1:
#    0.2686 pmol NPs ~ 370 pmol atoms
# DH2:
#    0.2696 pmol NPs ~ 371 pmol atoms
# DH3:
#    0.2712 pmol NPs ~ 374 pmol atoms
# ANH 12 [14409] raw:
#   First half: unrastered. 0.06545 pmol ~ 4.866 %
#   Second half: rastered.  0.1008 pmol ~ 1.873 %
#   Total coverage: 6.739 %
#   Total molage: 0.0907 pmol NPs ~ 500 pmol atoms
# ANH 13 [14415] : 0.1524 pmol NPs = 840 pmol atoms ~ 2.832 % --> 840/4 = 210 pmol atoms on sample
# ANH 14 [14422] : 0.2705 pmol NPs = 1491 pmol atoms ~ 5.028 % --> 1491/4 = 373 pmol atoms on sample
#                                                   /  20.11 %
