import LightwaveExplorer as Le
import numpy as np
import matplotlib.pyplot as plt

# help(LightwaveExplorer)


results_zip = "Data_BBO_sim/BBO_10percent_test.zip"

# results_zip = "Data_BBO/BBO_cont"

sim = Le.lightwaveExplorerResult(results_zip, loadFieldArray=True)

# print(vars(sim_30_microns).keys())   #we are interested in Ext_x, Ext_y, spectrum_x, spectrum_y

input_polarisation = sim.polarizationAngle1

print("pol input angle:" +str(input_polarisation*180/np.pi))


Nfreq = sim.Nfreq

print("Nfreq = " +str(Nfreq))

fStep = sim.fStep

print("fStep in THz = "+str(fStep*1e-12))


freq = np.linspace(0, Nfreq*fStep*1e-12, Nfreq)


beamwaist1 = sim.beamwaist1

print("input beam waist in microns of: "+str(beamwaist1*1e6))

crystalThickness = sim.crystalThickness

print("crystal thickness in microns: "+str(crystalThickness*1e6))







#####################for E_pol field in time - gives pulse shape

# print(type(sim_30_microns.Ext_x))   #(664,400) shape, 664 is for the time, 400 is for the x component , is a numpy array

E_x = sim.Ext_x  #so this is the total pulse, we just want to look at the on axis component!!
E_x_onaxis = E_x[:,200]


N_time = sim.Ntime
time_span = sim.timeSpan
time_step = sim.timeStep

print(N_time, time_span,time_step)

time_array = np.linspace(0,time_span*1e13, N_time)

# # plt.imshow(E_x)  #for full picture
plt.plot(time_array, E_x_onaxis)

plt.show()






###TO DO:  test SFG if pulse shape matters!! i.e. we can get more pulse overlap??

###################### for freq spectrum, try figure out if this correlates with the energy?


test = sim.spectrum_x* 1e12  # convert to J/THz

print(test.shape)

plt.plot(freq,test)
plt.ylabel("energy density / J per THz")
plt.xlabel("frequency / THz")

plt.show()


area = np.trapezoid(test, freq)

print("energy of x pol pulse in mJ: "+str(area*1e3))   # seems like is slightly inaccurate, but vaguely good enough for rough unit of measurement


test_y = sim.spectrum_y* 1e12

area_2 = np.trapezoid(test_y,freq)

print("energy of y pol pulse in mJ: " +str(area_2*1e3))



print("to check total energy in mJ: "+str(area*1e3 + area_2*1e3))  #seems like can be more than the total input - but we are just doing rough estimates...so is all fine, dont need to be super precise here!













