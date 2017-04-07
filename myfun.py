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
def run_script(IDENTIFIER, STRING):
    # Get parameters --------------------------------------
    SENSITIVITY_FILTER, SENSITIVITY_LIMIT, FIRST_LIMIT, PLOTTING, DEBUGGING, RADIUS_PARTICLE, RADIUS_APERTURE = get_parameters(STRING)
    # Load raw data and apply filters ---------------------
    time, current = get_filtered_data(IDENTIFIER, SENSITIVITY_FILTER, PLOTTING, DEBUGGING)
    # Analyze start ---------------------------------------
    LIMIT = FIRST_LIMIT
    first_region = get_first_region(time, current, LIMIT, SENSITIVITY_LIMIT, PLOTTING, DEBUGGING)
    # Extrapolate leak current measurements ---------------
    extrapolated_leak, deposition, leak = extrapolate_leak(time, current, first_region, LIMIT, SENSITIVITY_LIMIT, DEBUGGING)
    # Calaculate coverage and estimated remaining time ----
    get_info(time, current, deposition, extrapolated_leak, RADIUS_PARTICLE, RADIUS_APERTURE)
    # Plots -----------------------------------------------
    if PLOTTING and not DEBUGGING:
        plt.figure(1)
        for i in np.arange(deposition.viewkeys().__len__()):
            plt.plot(time[deposition[i]]/60, current[deposition[i]]*1e12, 'bo')
            plt.plot(time[deposition[i]]/60, extrapolated_leak[i]*1e12, 'yo')
        for i in np.arange(leak.viewkeys().__len__()):
            plt.plot(time[leak[i]]/60., current[leak[i]]*1e12, 'go')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Deposition current (pA)')
        plt.show()
    elif DEBUGGING:
        plt.show()
###########################################################
def get_parameters(STRING):
    STRING = STRING.split(';')
    #print('Analyzing deposition current using following parameters:\n')
    for i in STRING:
        x, y = i.split('=')
        if x == 'SENSITIVITY_FILTER':
            SENSITIVITY_FILTER = float(y)
            #print('\tSensitivity factor for noise filter: {}'.format(SENSITIVITY_FILTER))
        elif x == 'SENSITIVITY_LIMIT':
            SENSITIVITY_LIMIT = float(y)
            #print('\tSensitivity factor for leak detection: {}'.format(SENSITIVITY_LIMIT))
        elif x == 'FIRST_LIMIT':
            FIRST_LIMIT = float(y)*1e-12
            #print('\tMagnitude of first current: {} pA'.format(round(FIRST_LIMIT*1e12,2)))
        elif x == 'PLOT':
            PLOTTING = parseBoolString(y)
        elif x == 'PARTICLE_DIAMETER':
            RADIUS_PARTICLE = float(y)/2.*1e-9
            #print('\tDiameter of particles: {} nm'.format(RADIUS_PARTICLE*2e9))
        elif x == 'APERTURE_DIAMETER':
            RADIUS_APERTURE = float(y)/2.*1e-3
            #print('\tDiameter of aperture: {} mm'.format(RADIUS_APERTURE*2e3))
        elif x == 'DEBUG':
            DEBUGGING = parseBoolString(y)
    #if PLOTTING and not DEBUGGING:
    #    print('\tPlotting mode: normal')
    #elif DEBUGGING:
    #    print('\tPlotting mode: debugging/tuning')
    #else:
    #    print('\tPlotting mode: OFF')
    #print('')
    return SENSITIVITY_FILTER, SENSITIVITY_LIMIT, FIRST_LIMIT, PLOTTING, DEBUGGING, RADIUS_PARTICLE, RADIUS_APERTURE
#----------------------------------------------------------
def parseBoolString(STRING): # Interpret string as true or false
    return (STRING[0].upper()=='T') or (STRING[0].upper()=='ON') or (STRING[0].upper()=='Y')
###########################################################
def get_filtered_data(IDENTIFIER, SENSITIVITY_FILTER, PLOTTING, DEBUGGING):
    # Load raw data
    rtime, rcurrent = get_raw_data(IDENTIFIER)
    if PLOTTING and not DEBUGGING:
        plt.figure(1)
        plt.plot(rtime/60, rcurrent*1e12, 'ro-')
    # Compute averaged gradient
    averagedGradient = get_averagedGradient(rcurrent, width=2)
    std = get_gradient(rtime, rcurrent, DEBUGGING)
    # Apply noise filter
    #averagedGradient = smooth(abs(averagedGradient), 1)
    FILTER = filter_noise(smooth(abs(averagedGradient), 1), std, SENSITIVITY_FILTER)
    time = rtime[FILTER]
    current = rcurrent[FILTER]
    if DEBUGGING:
        plt.figure(1)
        plt.title('DEBUGGING: SENSITIVITY_FILTER')
        plt.plot(rtime/60, abs(averagedGradient)/std,'ro-')
        plt.plot(rtime/60, np.ones(len(rtime)), 'b-')
        plt.plot(rtime/60, np.ones(len(rtime))*SENSITIVITY_FILTER, 'b', linestyle='dashed')
    return time, current
#----------------------------------------------------------
def get_raw_data(IDENTIFIER):
    db = Cinfdata('omicron', use_caching=False, log_level='DEBUG')
    middata = db.get_data(IDENTIFIER)
    metadata = db.get_metadata(IDENTIFIER)
    print('\nLoading data from: {}'.format(metadata['Comment']))
    time = middata[:,0]
    current = middata[:,1]
    FILTER = filter_overflow(current)
    return time[FILTER], current[FILTER]
#----------------------------------------------------------
def filter_overflow(current): # Remove overflowed data where I > 0
    return np.where(current < 1)
#----------------------------------------------------------
def get_averagedGradient(current, width=2):
    ### INTERPRETATION: ###
    # gradient[i] related to time[i] is the change at time[i] to time[i+1]
    # gradient[i] related to time[i+1] is the change at time[i+1] from time[i]
    # averagedGradient is the gradient averaged over *width* neighbors
    gradient = current[1:] - current[0:-1]
    NUM = len(current)
    averagedGradient = np.zeros(NUM)
    # calculate middle
    DENOMINATOR = 2*width + 1
    for i in range(DENOMINATOR-1):
        averagedGradient[width:-width-1] += gradient[i:i-2*width]/DENOMINATOR
    averagedGradient[width:-width-1] += gradient[i+1:]/DENOMINATOR
    # append ends
    averagedGradient[0:width] = averagedGradient[width]
    averagedGradient[-width-1:] = averagedGradient[-width-2]
    return averagedGradient
#----------------------------------------------------------
def get_gradient(time, current, DEBUGGING):
    counter = 0
    OFFSET = 15
    while True:
        gradient = abs(current[1:OFFSET+counter]-current[0:OFFSET+counter-1])
        stdGradient = gradient.std()
        averagedGradient = get_averagedGradient(current[0:OFFSET+counter+6])
        counter += 1
        if abs(averagedGradient[OFFSET+counter]) > stdGradient:
            print('Gradient computed up to index {}'.format(OFFSET+counter))
            break
    if DEBUGGING:
        counter = 0
        while True:
            gradient = current[1:OFFSET+counter]-current[0:OFFSET+counter-1]
            new_stdGradient = gradient.std()
            averagedGradient = get_averagedGradient(current[0:OFFSET+counter+6])
            plt.figure(1)
            plt.plot(time[counter+OFFSET]/60, new_stdGradient/stdGradient, 'bo')
            counter += 1
            if abs(averagedGradient[OFFSET+counter]) > stdGradient:
                break
    return stdGradient
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
def filter_noise(reducedGradient, stdGradient, SENSITIVITY_FILTER): # Filter noise
    return np.where(abs(reducedGradient) < stdGradient*SENSITIVITY_FILTER)
###########################################################
def get_first_region(time, data, FIRST_LIMIT, SENSITIVITY_LIMIT, PLOTTING, DEBUGGING):
    counter = 0
    OFFSET = 15
    try:
        while True:
            index = OFFSET + counter
            if data[index+1] - data[index] > FIRST_LIMIT*SENSITIVITY_LIMIT:
                first_region = 'deposition'
                break
            elif data[index+1] - data[index] < -FIRST_LIMIT*SENSITIVITY_LIMIT:
                first_region = 'background'
                break
            counter += 1
        print('First region detected: "{0}" at filtered index {1} = {2} minutes'.format(first_region, index, round(time[index]/60, 2)))
        plot_FIRST_LIMIT(time, data, 0, index, first_region, FIRST_LIMIT, SENSITIVITY_LIMIT, PLOTTING, DEBUGGING)
    except IndexError:
        print('   *** INDEXERROR: try changing FIRST_LIMIT or SENSITIVITY_LIMIT   *** ')
        plot_FIRST_LIMIT(time, data, 0, -2, 'background', FIRST_LIMIT, SENSITIVITY_LIMIT, PLOTTING, DEBUGGING)
        plot_FIRST_LIMIT(time, data, 0, -2, 'deposition', FIRST_LIMIT, SENSITIVITY_LIMIT, PLOTTING, DEBUGGING)
        plt.show()
    return first_region
#----------------------------------------------------------
def plot_FIRST_LIMIT(time, data, index, index_X, first_region, FIRST_LIMIT, SENSITIVITY_LIMIT, PLOTTING, DEBUGGING):
    if PLOTTING and not DEBUGGING:
        plt.figure(1)
        LIMIT = FIRST_LIMIT*1e12
        x = np.array([0, time[index_X+1], time[index_X+1], 0])/60
        if first_region == 'deposition':
            y = np.array([0,0,LIMIT,LIMIT])+data[index]*1e12
            plt.fill(x, y, color='k', alpha=0.2)
            LIMIT = LIMIT*SENSITIVITY_LIMIT
            y = np.array([0,0,LIMIT,LIMIT])+data[index]*1e12
            plt.fill(x, y, color='k', alpha=0.5)
        elif first_region == 'background':
            y = np.array([0,0,-LIMIT,-LIMIT])+data[index]*1e12
            plt.fill(x, y, color='k', alpha=0.2)
            LIMIT = LIMIT*SENSITIVITY_LIMIT
            y = np.array([0,0,-LIMIT,-LIMIT])+data[index]*1e12
            plt.fill(x, y, color='k', alpha=0.5)
###########################################################
def extrapolate_leak(time, data, first_region, LIMIT, SENSITIVITY_LIMIT, DEBUGGING):
    deposition, leak, last_region = get_regions(time, data, first_region, LIMIT, SENSITIVITY_LIMIT, DEBUGGING)
    counter_current = deposition.viewkeys().__len__()
    extrapolated_leak = dict()
    if first_region == 'deposition':
        for i in np.arange(counter_current):
            if i == 0:
                a, b = 0, data[leak[i]].mean()
            elif (i == counter_current - 1) and (last_region == 'deposition'):
                a, b = 0, data[leak[i-1]].mean()
            else:
                a = (data[leak[i]].mean() - data[leak[i-1]].mean())/(time[leak[i][-1]]-time[leak[i-1][0]])
                b = (data[leak[i]].mean() - a*time[leak[i][-1]])
            extrapolated_leak[i] = a*time[deposition[i]] + b
    else:
        for i in np.arange(counter_current):
            if (i == counter_current - 1) and (last_region == 'deposition'):
                a, b = 0, data[leak[i]].mean()
            else:
                a = (data[leak[i]].mean() - data[leak[i+1]].mean())/(time[leak[i][-1]]-time[leak[i+1][0]])
                b = (data[leak[i]].mean() - a*time[leak[i][-1]])
            extrapolated_leak[i] = a*time[deposition[i]] + b
    return extrapolated_leak, deposition, leak
#----------------------------------------------------------
def get_regions(time, data, first_region, LIMIT, SENSITIVITY_LIMIT, DEBUGGING):
    # Automatically separate deposition current into regions
    # of background and actual deposition current
    deposition, leak = dict(), dict()
    counter_current, counter_leak = 0, 0
    if first_region == 'background':
        depo = False
    else:
        depo = True
    previous = -1
    if DEBUGGING:
        plt.figure(2)
        plt.title('DEBUGGING: SENSITIVITY_LIMIT')
        plt.plot(time[1:]/60, abs(data[1:]-data[0:-1]), 'bo-')
        plt.plot(time/60, np.ones(len(time))*LIMIT, 'r-')
    # Index i represents last index of region
    for i in np.arange(len(data)-1):
        # Detect end of leak measurement
        if ((data[i+1] < data[i] - SENSITIVITY_LIMIT*LIMIT) and not depo):
            leak[counter_leak] = np.arange(previous+1,i+1)
            previous = i
            counter_leak += 1
            depo = True
            if DEBUGGING and (counter_leak == 1):
                plt.plot(time[i:]/60, np.ones(len(time[i:]))*LIMIT*SENSITIVITY_LIMIT, 'r', linestyle='dashed')
        # Detect end of current measurement
        elif ((data[i+1] > data[i] + SENSITIVITY_LIMIT*LIMIT) and depo):
            deposition[counter_current] = np.arange(previous+1,i+1)
            previous = i
            counter_current += 1
            depo = False
            LIMIT = renew_LIMIT(data, i)
            if DEBUGGING:
                plt.plot(time[i:]/60, np.ones(len(time[i:]))*SENSITIVITY_LIMIT*LIMIT, 'r', linestyle='dashed')
    # End of measurement
    if depo:
        last_region = 'deposition'
        deposition[counter_current] = np.arange(previous+1,len(data))
    else:
        last_region = 'background'
        leak[counter_leak] = np.arange(previous+1,len(data))
    print('\nRegions detected:\n\t{} deposition currents and\n\t{} leak currents\n'.format(deposition.viewkeys().__len__(), leak.viewkeys().__len__()))
    return deposition, leak, last_region
#----------------------------------------------------------
def renew_LIMIT(data, index):
    left = data[index-5:index+1]
    left = sum(left)/len(left)
    right = data[index+1:index+7]
    right = sum(right)/len(right)
    return abs(left-right)
###########################################################
def get_info(time, current, deposition, extrapolated_leak, RADIUS_PARTICLE, RADIUS_APERTURE):
    number_of_charges, charges_per_second = integrate_current(time, current, deposition, extrapolated_leak)
    TARGET, coverage, time_remaining_minutes, time_remaining_seconds = calculate_coverage(number_of_charges, charges_per_second, RADIUS_PARTICLE, RADIUS_APERTURE)
    print('Calculated coverage:\n   ***   {0} percent   ***\n'.format(round(coverage,3)))
    print('Estimated remaining time to {0} percent coverage:\n   ***   {1} minutes and {2} seconds   ***'.format(TARGET, time_remaining_minutes, round(time_remaining_seconds, 1)))
#----------------------------------------------------------
def integrate_current(time, current, deposition, extrapolated_leak):
    # Return number of nanoparticles
    NUM = deposition.viewkeys().__len__()
    integration = 0
    for i in np.arange(NUM):
        actual_current = extrapolated_leak[i] - current[deposition[i]]
        actual_current = (actual_current[0:-1] + actual_current[1:])/2.
        integration += sum(actual_current*np.diff(time[deposition[i]]))
    e = 1.602e-19
    print('Number of charges presently measured: {0:.4} pmol'.format(integration/e/6.022e23*1e12))
    #print('Charges per second: {}'.format(actual_current.mean()/e))
    return integration/e, actual_current.mean()/e
#----------------------------------------------------------
def calculate_coverage(number_of_charges, charges_per_second, RADIUS_PARTICLE, RADIUS_APERTURE):
    area_particle = np.pi * ((RADIUS_PARTICLE)**2)
    projected_area = area_particle * number_of_charges
    area_aperture = np.pi * ((RADIUS_APERTURE)**2)
    coverage = projected_area/area_aperture * 100
    for TARGET in [5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300]:
        if coverage < TARGET:
            total_number_of_charges = TARGET/100.*area_aperture/area_particle
            relative = charges_per_second*0.128
            #print('Single step percent of current measurement: {0:.4}'.format(relative/number_of_charges*100))
            #print('Single step percent of target measurement: {0:.4}'.format(relative/total_number_of_charges*100))
            time_remaining = (total_number_of_charges - number_of_charges)/charges_per_second
            time_remaining_minutes = int(time_remaining)/60
            time_remaining_seconds = time_remaining - time_remaining_minutes*60
            break
    return TARGET, coverage, time_remaining_minutes, time_remaining_seconds
################################################################
### MAIN ###
################################################################
if __name__ == '__main__':
    try:
        run_script(11604, 'PARTICLE_DIAMETER=3.5;APERTURE_DIAMETER=9.;FIRST_LIMIT=1.0;SENSITIVITY_LIMIT=.6;SENSITIVITY_FILTER=1.;PLOT=True;DEBUG=rue')
    except:
        print('***\nSomething is wrong: Check input parameters or try debugging mode!!\n***')
        plt.show()
        raise
    print('--- END OF SCRIPT ---')
