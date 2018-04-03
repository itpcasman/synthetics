# -*- coding: utf-8 -*-

"""
Create synthetic seismics from petrophysical data (MSCL) of piston cores.
Parameters used: bulk density (rho), slowness (1/p-wave velocity).
Note: The code needs to be readjusted if the synthetic seismogrsm of a core not starting at the sea floor needs to be calculated.
Requires python >= 3.5

@author: Andreas Paul
@year: 2015 - 2018
"""

#
# Import modules required
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#
# Function definitions
#

def w_pressure(depth): 
    """ Calculate water pressure at different depths in kPa """
    w_pres = 1.03e3 * 9.8 * depth
    w_pres = w_pres + 1.01e5
    w_pres = w_pres / 1000
    return w_pres

def w_acoustic_vel(T,S,Z,lat):
    """ Calculate acoustic velocity of water dependent on water depth, temperature, salinity and latitude. After Leroy et al. (2008) J. Acoust. Soc. Am. 124(5). """
    w_ac_vel = 1402.5 + 5 * T - 5.44e-2 * T**2 + 2.1e-4 * T**3 + 1.33 * S - 1.23e-2 * S * T + 8.7e-5 * S * T**2 + 1.56e-2 * Z + 2.55e-7 * Z**2 - 7.3e-12 * Z**3 + 1.2e-6 * Z * (lat - 45) - 9.5e-13 * T * Z**3 + 3e-7 * T**2 * Z + 1.43e-5 * S * Z    
    return w_ac_vel

def ricker(f, length, dt):
    """ Define Ricker wavelet """
    t = np.linspace(-length / 2, (length-dt) / 2, length / dt)
    y = (1. - 2.*(np.pi**2)*(f**2)*(t**2))*np.exp(-(np.pi**2)*(f**2)*(t**2))
    return t, y


#
# Data import
# Important: Before importing, the data has to be cleaned of 'bad' data such as rows with p-wave amplitudes < 95 or spikes in density related to core section transitions.
# Cleaning can be done in Excel or with Pandas but is not implemented/shown here. Edit for using your own .csv files.
# Data import
#
    
labels = np.genfromtxt('1150_petrophysical.csv', delimiter=';', dtype=str)[0,:]
petrophys_1150 = np.genfromtxt('1150_petrophysical.csv', delimiter=";", skip_header=1)


#
# Extract and define variables from the variable petrophys_1150; calculate tvdss (z) and dts. Edit for your own file.
#

td = petrophys_1150[:,0]                        # Total depth below sea-floor [m]
ap = petrophys_1150[:,5]                        # P-wave amplitude [0-100]
vp = petrophys_1150[:,6]                        # P-wave velocity [m/s]
rho = petrophys_1150[:,7]                       # Bulk density [g/cm^3]
phi = petrophys_1150[:,8]                       # Fractional porosity [0-1]
z = td[0 :] + 245                               # Total vertical depth below sea level [m]
dts = 1e6 / vp                                  # Convert p-wave velocity to slowness [us/m]


#
# Depth - Time conversion and velocity model
#

core_depth = 245    # Depth at coring location in meters 
water_pres = w_pressure(core_depth)# Calculate water pressure at core depth
v_wat_core_depth = w_acoustic_vel(20,3.5,core_depth,4)  # Water acoustic vel at core depth (20 C, 3.5% salinity, core depth, 4 deg lat)
v_wat_surf = w_acoustic_vel(28,3.5,0,4) # Water acoustic velocity at surface

# Calculate TWT in water in s (sea-level --> sea floor --> sea-level)
owt_water = np.abs(core_depth) / ((v_wat_core_depth + v_wat_surf) / 2)
twt_water = 2.0 * (np.abs(core_depth) / ((v_wat_core_depth + v_wat_surf) / 2))
log_start_time = owt_water 
log_start_time_ms = owt_water * 1000
log_start_time_twt = twt_water

# Convert to int (round)
v_wat_core_depth = int(v_wat_core_depth)
v_wat_surf = int(v_wat_surf)


#
# Clean and despike the data, remove extreme values, interpolate gaps (Not used here)
#

rho_desp = rho
dts_desp = dts


#
# Conversion and interpolation (not used!)
#

# td_int = np.arange(0, 14.50, 0.001)
# dts_int = np.interp(td_int, td, dts_desp)
# rho_int = np.interp(td_int, td, rho_desp)
# z_int = td_int[0 :] + 245


#
# Two-way-time to depth relationship
#

scaled_dts = 0.01 * dts_desp / 1e6
tcum = 2 * np.cumsum(scaled_dts)
tdr = (tcum + log_start_time) * 2


#
# Acoustic impedance
#

imp = dts_desp * rho_desp


#
# Reflection coefficients
#

rc = (imp[1:] - imp[:-1]) / (imp[1:] + imp[:-1])


#
# Create the synthetics
#

f = 18  # Frequency
tw, w = ricker(f, length = 0.512, dt = 0.001)
synth = np.convolve(w, rc, mode='same')


#
# Plotting vs depth, time and seismic
#

f = plt.figure(figsize=[20,10],facecolor='white')
gs = gridspec.GridSpec(1, 7, width_ratios = [1.25, 1.25, 0.5, 1.25, 0.75, 0.75, 8]) # The grid created here originally included a plot of the actual seismic section which was removed here. One can adjust the gridspec to only have 5 columns(cells). The width_ratios need to be adjusted accordingly.


#
# depth domain
#

axa = plt.subplot(gs[0])
axa.plot(dts, td,'k', alpha=0.75, lw=0.25)
axa.set_title('Travel Time', fontsize=10)
axa.set_ylabel('Depth (m) ', fontsize = '10' )
axa.set_xlabel('P-wave slowness 'r'$[\mu s/m]$', fontsize = '10')
axa.set_ylim(0, 14.5)
# axa.set_xticklabels('')
axa.tick_params(labelsize=8)
plt.setp(axa.xaxis.get_majorticklabels(), rotation=45)
axa.invert_yaxis()
axa.grid()

gs.update(wspace=0.1)

axb = plt.subplot(gs[1])
axb.plot(tdr, td, 'k', alpha = 0.75)
axb.set_title('Interval travel time', fontsize=10)
axb.set_xlabel('TWT ' + '$[s]$', fontsize = '10')
axb.set_ylim(0, 14.5)
axb.invert_yaxis()
axb.tick_params(labelsize=8)
axb.set_yticklabels('')
plt.setp(axb.xaxis.get_majorticklabels(), rotation=45)
# axb.set_xticklabels('')

axb.grid()

#
# time domain
#

# white space between depth and time plots
axoff = plt.subplot(gs[2])
axoff.set_axis_off()

axc = plt.subplot(gs[3])
axc.plot(imp, tdr, 'k', alpha=0.75)
axc.set_title('Impedance', fontsize=10)
axc.set_ylabel('TWT '+ '$[s]$', fontsize = '10' )
axc.set_xlabel(r'$kg / m^2s^2$ ', fontsize = '10')
axc.tick_params(labelsize=8)
plt.setp(axc.xaxis.get_majorticklabels(), rotation=45)
# axc.set_ylim(0.32, 0.72)
axc.set_ylim(0.32, 0.36)
axc.set_xlim(1200, 1400)
axc.invert_yaxis()
axc.grid()

axd = plt.subplot(gs[4])
axd.hlines(tdr[:-1], 0, rc, color='k', lw = 1)                    # Stems
axd.plot([0, 0], [tdr.min(), tdr.max()], '-', c='k', alpha = 0.5) # Middle bar
axd.set_title('Reflectivity', fontsize=10)
axd.set_xlabel('', fontsize = '12')
# axd.set_ylim(0.32, 0.72)
axd.set_ylim(0.32,0.36)
axd.set_xlim(-0.025, 0.025)
axd.set_yticklabels('')
axd.set_axis_off()
axd.invert_yaxis()

axe = plt.subplot(gs[5])
axe.plot( synth, tdr[:-1], color='none')
axe.set_title('Synthetics', fontsize=10)
axe.fill_betweenx(tdr[:-1], synth, 0, synth > 0.0, color='k', alpha = 1.0)
axe.fill_betweenx(tdr[:-1], synth, 0, synth < 0.0, color='red', alpha = 1.0)
axe.set_xlabel('', fontsize = '12')
# axe.set_ylim(0.32, 0.72)
axe.set_ylim(0.32,0.36)
axe.set_xlim(-0.05, 0.05)
axe.set_yticklabels('')
axe.set_axis_off()
axe.invert_yaxis()

plt.show()


