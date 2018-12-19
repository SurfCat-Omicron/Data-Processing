import sys
import matplotlib.pyplot as plt
import numpy as np

"""
Basic version copied over from "ISS_resolution.py" which gives the plot used in my master thesis.

Since then, mass 16, 18, and 19 has been added and energies changed.
"""

###################
### Definitions ###
###################
# Plot relationship between kinetic energy and target mass
def plot_relationship(E0=1000):
    fig = plt.figure()
    global ax
    ax = fig.add_subplot(111)
    # Plot helium
    mass_range_He = np.arange(3,300,0.1)
    ax.plot(mass_range_He, mass2energy(mass_range_He, mi=3, E0=E0), 'k-')
    #ax.text(15, 1300, 'Helium-3', verticalalignment='top', horizontalalignment='center')
    # Plot helium
    mass_range_He = np.arange(4,300,0.1)
    ax.plot(mass_range_He, mass2energy(mass_range_He, mi=4, E0=E0), 'k-')
    #ax.text(15, 400, 'Helium-4', verticalalignment='top', horizontalalignment='center')
    # Plot Neon
    mass_range_Ne = np.arange(20,300,5)
    ax.plot(mass_range_Ne, mass2energy(mass_range_Ne, mi=20, E0=E0), 'k-.')
    ax.text(26, 150, 'Neon', verticalalignment='top', horizontalalignment='center')
    # Plot Argon
    mass_range_Ar = np.arange(40,300,5)
    ax.plot(mass_range_Ar, mass2energy(mass_range_Ar, mi=40, E0=E0), 'k:')
    ax.text(70, 80, 'Argon', verticalalignment='top', horizontalalignment='left')
    # Axis labels
    ax.set_xlabel('Target mass (amu)')
    ax.set_ylabel('Kinetic energy (eV)')
    [x1, x2, y1, y2] = ax.axis()
    alpha=0.4
    # Add 3d-metals
    xmin, xmax = 45, 65
    ax.fill([xmin,xmin,0,0,xmax,xmax],
             mass2energy(np.array([4,xmin,xmin,xmax,xmax,4]), E0=E0),
             color='k', alpha=alpha)
    # Add masses 16, 18, 19
    for m in [3, 4]:
        for m2 in [16, 18, 19]:
            ax.plot(m2, mass2energy(m2, mi=m, E0=E0), 'ko', markersize=8)
    #ax.text((xmax+xmin)/2., 300,
    #        '3d metals',
    #        fontweight='bold',
    #        horizontalalignment='center')
    #ax.annotate('3d metals',
    #            xy=((xmax+xmin)/2., 400),
    #            xytext=((xmax+2*xmin)/3., 300),
    #            arrowprops=dict(arrowstyle='->'),
    #            color='k',
    #            fontweight='bold',
    #            horizontalalignment='right')
    # Add 4d-metals
    xmin, xmax = 89, 112
    ax.fill([xmin,xmin,0,0,xmax,xmax],
             mass2energy(np.array([4,xmin,xmin,xmax,xmax,4]), E0=E0),
             color='k', alpha=alpha)
    #ax.text((xmax+xmin)/2., 400,
    #        '4d metals',
    #        fontweight='bold',
    #        horizontalalignment='center')
    #ax.annotate('4d metals',
    #            xy=((xmax+xmin)/2., 500),
    #            xytext=((xmax+2*xmin)/3., 400),
    #            arrowprops=dict(arrowstyle='->'),
    #            color='k',
    #            fontweight='bold',
    #            horizontalalignment='right')
    # Add 5d-metals
    xmin, xmax = 178, 201
    ax.fill([xmin,xmin,0,0,xmax,xmax],
             mass2energy(np.array([4,xmin,xmin,xmax,xmax,4]), E0=E0),
             color='k', alpha=alpha)
    #ax.text((xmax+xmin)/2., 600,
    #        '5d metals',
    #        fontweight='bold',
    #        horizontalalignment='center')
        #ax.annotate('5d metals',
    #            xy=((xmax+xmin)/2., 700),
    #            xytext=((2*xmax+xmin)/3., 600),
    #            arrowprops=dict(arrowstyle='->'),
    #            color='k',
    #            fontweight='bold',
    #            horizontalalignment='left')
    # Add Lanthanides
    xmin, xmax = 139, 175
    ax.fill([xmin,xmin,0,0,xmax,xmax],
             mass2energy(np.array([4,xmin,xmin,xmax,xmax,4]), E0=E0),
             color='y', alpha=alpha)
    #ax.text((xmax+xmin)/2., 500,
    #        'Lanthanides',
    #        fontweight='bold',
    #        horizontalalignment='center')
    #ax.annotate('Lanthanides',
    #            xy=((xmax+xmin)/2., 600),
    #            xytext=((xmax+2*xmin)/3., 500),
    #            arrowprops=dict(arrowstyle='->'),
    #            color='k',
    #            fontweight='bold',
    #            horizontalalignment='right')
    # Add Actinides
    xmin, xmax = 227, 262
    ax.fill([xmin,xmin,0,0,xmax,xmax],
             mass2energy(np.array([4,xmin,xmin,xmax,xmax,4]), E0=E0),
             color='y', alpha=alpha)
    #ax.text((xmax+xmin)/2., 700,
    #        'Actinides',
    #        fontweight='bold',
    #        horizontalalignment='center')
    #ax.annotate('Actinides',
    #            xy=((xmax+xmin)/2., 800),
    #            xytext=((xmax+2*xmin)/3., 700),
    #            arrowprops=dict(arrowstyle='->'),
    #            color='k',
    #            fontweight='bold',
    #            horizontalalignment='right')
    #draw_infobox(ax, orientation='lower right')
    #ax.set_yticks([0,200,400,600,800,1000, 1200, 1400, 1600, 1800, 2000])
    #plt.show()

# Convert a mass to corresponding reflected energy
def mass2energy(ms, E0=1000, mi=4, theta=146.7):
    theta = theta * np.pi/180
    return E0 * ( (mi*np.cos(theta) + np.sqrt(ms**2 - mi**2*np.sin(theta)**2))/(ms + mi) )**2

# Draw settings of the ISS in upper corner
def draw_infobox(axis, orientation='upper right'):
    if orientation == 'upper right':
        x, y = 0.8, 0.95
    elif orientation == 'upper left':
        x, y = 0.05, 0.95
    elif orientation == 'lower right':
        x, y = 0.8, 0.3
    elif orientation == 'lower left':
        x, y = 0.05, 0.3
    axis.text(x, y, '$Incident$ $ions$:\n\t$\\theta = 146.7^\circ$\n\t$Type:$ $He^+$\n\t$E_0 = 1000$ $eV$', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes)

# Make two axes with energies and reflected mass
def convert_energies(axis, elements):
    ax2=axis.twiny()
    ax2.set_xlabel('Kinetic energy (eV)')
    ax2.set_xlim(axis.get_xlim())
    ax2.set_xticks(axis.get_xticks())
    ax2.set_xticklabels([str(x*200) for x in range(6)])
    axis.set_xlabel('Target mass (amu)')
    axis.set_xticks([mass2energy(m) for m in elements])
    axis.set_xticklabels(elements)

# Find start points
def find_start(contents):
    return [i for i in range(len(contents)) if contents[i] == 'Energy\tCounts\r\n']

# Open files to get data
def open_file(folder, files):
    f = open(folder+'/'+files)
    data = f.readlines()
    f.close()
    return data

# Retrieve dwell times
def get_dwells(data, i):
    line = data[i-4].split('\t')
    for j in range(len(line)):
        if line[j].lower() == 'dwell':
            return float(data[i-3].split('\t')[j])

# Copy data
def get_data(data, end, i):
    energy, counts = np.zeros(end)+1., np.zeros(end)+1.
    counter = 0
    for j in np.arange(end) + i + 1:
        middata = data[j].rstrip('\r\n').split('\t')
        energy[counter], counts[counter] = float(middata[0]), float(middata[1])
        counter += 1
    return energy, counts

# Plot element markers
def draw_guide(axis, element, mass, scale=.57):
    x = mass2energy(mass)
    [x1, x2, y1, y2] = axis.axis()
    axis.axvline(x=x, ymin=0, ymax=scale-0.02, color='k', linestyle='dotted')
    axis.text(x, scale*y2, element, horizontalalignment='left')

# Plot all consecutive measurements
def plotall(axis, start, energy, counts, dwell, color):
    n = len(start)
    alpha = (np.arange(0,n)+1)/float(n*2)
    alpha[-1] = 1.0
    counter = 0
    for run in start:
        axis.plot(energy[run], counts[run]/dwell[run], marker='o', linestyle='', color=color, alpha=alpha[counter])
        counter += 1

# Plot single run (measurement)
def plot1(axis, run, start, energy, counts, dwell, color, alpha=1.0):
    key = start[run]
    axis.plot(energy[key], counts[key]/dwell[key], marker='o', linestyle='', color=color, alpha=alpha)

# Auto-annotate: Add peak data desciptors to maximum of peak
# Note: cps = counts/dwell
def autoannotate(axis, key, run, energy, cps, color, horizontalalignment='center'):
    pointer = cps.argmax()
    x = energy[pointer]
    y = cps[pointer]
    [x1, x2, y1, y2] = axis.axis()
    ynew = y+y2*0.06
    if ynew > y2:
        axis.axis([x1, x2, y1, ynew*1.05])
    axis.annotate(key,
                  xy=(x-x2*0.01, y),
                  xytext=(x-x2*0.10, ynew),
                  arrowprops=dict(arrowstyle='->'),
                  color=color,
                  fontweight='bold',
                  horizontalalignment=horizontalalignment)
#autoannotate(ax, key, run, energy[key], counts[key][run]/dwell[key][run], color[key],'left')

# Create an enlarged section of the data in a boxed plot
def zoomplot(fig, energy, cps, color, pos=[0.25,0.4,0.3,0.3], zoom=[600,1000,0,500]):
    ax = fig.add_subplot(211)
    ax.set_position(pos = pos)
    ax.plot(energy, cps, color=color)
    ax.axis(zoom)
    return ax

# Create an enlarged section of the data
def zoom(ax, energy, cps, color, alpha=0.5, zoom=600, factor=100, axis=True):
    [x1, x2, y1, y2] = ax.axis()
    n = np.where(energy==round(zoom))[0][0]
    ax.plot(energy[n:], cps[n:]*float(factor), '-', color=color, alpha=1)
    height = (round(max(cps[n:])*factor/1000)+1)*1000
    if axis:
        ax.text(energy[n], (round(max(cps[n:n+60])*factor/1000)+1)*1000, 'x'+str(factor), color='k', horizontalalignment='right')
        ax.arrow(energy[n]-50, 0, 0, height, head_width=0.01*(x2-x1), head_length=0.01*(y2-y1), fc='k', ec='k')
        ax.arrow(energy[n]-50-0.005*(x2-x1), (round(max(cps[n:])*factor/100))*100, 0.01*(x2-x1), 0, ec='k')
        ax.text(energy[n]-50-0.01*(x2-x1), (round(max(cps[n:])*factor/100))*100, str(int(round(max(cps[n:])*factor/100))*100/factor), horizontalalignment='right', verticalalignment='center', color='k')

##############
### Inputs ###
##############
for i in [1000, 2000, 3000]:
    #plt.figure
    plot_relationship(E0=i) # Uncomment to generate plot of the mass resolution of ISS with the given settings. Close the first plot without maximizing it first lest the infobox with ISS settings will be drawn in the wrong spot.
plt.show()
