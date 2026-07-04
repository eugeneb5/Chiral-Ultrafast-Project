import LightwaveExplorer as Le
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
import random

# help(Le.lightwaveExplorerResult)

def define_pulse(zipname, read_char = True, pol_y = False):

    sim = Le.lightwaveExplorerResult(zipname, loadFieldArray=True)

    
                
    input_polarisation = sim.polarizationAngle1


    Nfreq = sim.Nfreq
    
    fStep = sim.fStep

    # freq = np.linspace(0, Nfreq*fStep*1e-12, Nfreq)
    beamwaist1 = sim.beamwaist1

    crystalThickness = sim.crystalThickness

    if pol_y:

        E_y = sim.Ext_y

        midpoint = E_y.shape[1]//2
        E_onaxis = E_y[:,midpoint]
        full_E_field = E_y
    else:

        E_x = sim.Ext_x  #so this is the total pulse, we just want to look at the on axis component!!

        midpoint = E_x.shape[1]//2
        E_onaxis = E_x[:,midpoint]
        full_E_field = E_x

    N_time = sim.Ntime
    time_span = sim.timeSpan*1e15
    time_step = sim.timeStep  #is in s
    time_step_fs = time_step*1e15  #adjusted to units of fs!
    
    print(N_time, time_span,time_step)

    spatial_step = sim.spatialStep
    # spatial_height = sim.spatialHeight  #this is y var - so ignore since used 2D simulation!
    spatial_width = sim.spatialWidth*1e6  #this is x var - so use this!, is in microns too
    N_space = sim.Nspace

    # time_array = np.linspace(0,time_span, N_time)
    time_array = np.linspace(-time_span//2 , time_span//2 , N_time)
    space_array = np.linspace(-spatial_width//2, spatial_width//2, N_space)
    # print("spatial height and width: "+str(spatial_height)+str(spatial_width))

    if read_char:
        print("Nfreq = " +str(Nfreq))
        print("pol input angle:" +str(input_polarisation*180/np.pi))  #in degrees
        print("fStep in THz = "+str(fStep*1e-12))
        print("input beam waist in microns of: "+str(beamwaist1*1e6))
        print("crystal thickness in microns: "+str(crystalThickness*1e6))
        print("timestep in fs: "+str())
        
        plt.plot(time_array, E_onaxis)
        plt.show()


    return time_array, E_onaxis, time_step_fs, full_E_field, space_array





results_zip = "Data_BBO_sim/BBO_USE_THIS_TEST.zip"

pre_BBO_pulse = "Data_BBO_sim/BBO_vacuum_USE_THIS.zip"


time_array_BBO_pulse, E_onaxis_BBO_pulse, time_step_fs,full_E_field,space_array = define_pulse(results_zip)

time_array_pre_BBO_pulse, E_onaxis_pre_BBO_pulse,_,_,_ = define_pulse(pre_BBO_pulse, pol_y=True)







def waist_estimation(intensity_2D, time_array):

    fluence = np.trapezoid()

####try define envelope function?



def envelope_intensity_temporal_function(time_array, E_x_onaxis, n = 1, plot_mode = False):   #default n is for air, need to adjust for chiral solutions... we also are assuming linewidth is appropriately small so n variation negligible
    E_complex = hilbert(E_x_onaxis)

    Intensity_envelope = 0.5* 3e8 * n * 8.85e-12 * np.abs(E_complex)**2  #we keep as m^-2 here for consistency, change to cm^-2 by *1e4!

    if plot_mode:

        plt.plot(time_array,Intensity_envelope)


        plt.show()

    return Intensity_envelope

def intensity_2D(space_array,time_array,E_x, n = 1, plot_mode = False):   #default n is for air, need to adjust for chiral solutions... we also are assuming linewidth is appropriately small so n variation negligible
    E_complex = hilbert(E_x, axis = 0)

    Intensity_2D = 0.5* 3e8 * n * 8.85e-12 * np.abs(E_complex)**2  #we keep as m^-2 here for consistency, change to cm^-2 by *1e4!
    # dimensions = np.shape(Intensity_2D)
    # time_dim = dimensions[0]
    # space_dim = dimensions[1]
    plot_extent = [space_array[0],space_array[-1],time_array[0],time_array[-1]]
    if plot_mode:

        
        im = plt.imshow(Intensity_2D, extent=plot_extent, origin='lower', aspect='auto', cmap='viridis')


        plt.xlabel("x /microns")
        plt.ylabel("t /fs")


        cbar = plt.colorbar(im)
        cbar.set_label("Intensity")

        plt.show() 

    return Intensity_2D


Intensity_2D = intensity_2D(space_array, time_array_BBO_pulse,full_E_field,plot_mode = True)





I_1 = envelope_intensity_temporal_function(time_array_BBO_pulse,E_onaxis_BBO_pulse, plot_mode = False)

I_2 = envelope_intensity_temporal_function(time_array_pre_BBO_pulse,E_onaxis_pre_BBO_pulse, plot_mode=False)





def overlap_percentage(I_1, I_2):  #using cauchy schwarz?!

    overlap_pulse = I_1 * I_2 

    numerator = np.trapezoid(overlap_pulse)

    denominator = ((np.trapezoid(I_1**2))*(np.trapezoid(I_2**2)))**0.5

    pulse_overlap = numerator/denominator
    
    test_overlap = I_1*I_1

    test_numerator = np.trapezoid(test_overlap)

    test_denominator = ((np.trapezoid(I_1**2))*(np.trapezoid(I_1**2)))**0.5

    pulse_test = test_numerator/test_denominator

    return pulse_overlap/pulse_test   #pulse_test is 1 anyways, so no need!!




def time_gradation(delta_x):

    return delta_x*2/3e8




def generate_random_phase_delay(I_1,I_2,dt_fs,original_time_array,start=1,end=100, plot_mode = False):  #start, end parameters are in microns!!

    path_length_delay = time_gradation(random.randint(start,end)*1e-6)*1e15  #is in fs!
    # path_length_delay = -50
    # print("path length difference generated, in fs: "+str(path_length_delay))

    time_delay_array = original_time_array + path_length_delay

    array_increase = path_length_delay / dt_fs

    t_min = min(original_time_array[0], time_delay_array[0])
    t_max = max(original_time_array[-1], time_delay_array[-1])
    num_steps = int(np.round((t_max - t_min) / dt_fs)) + 1
    expanded_time_array = np.linspace(t_min, t_max, num_steps)
    I_1_expanded = np.interp(expanded_time_array, time_delay_array, I_1, left=0.0, right=0.0)   #need to understand how this works
    I_2_expanded = np.interp(expanded_time_array, original_time_array, I_2, left=0.0, right=0.0)
    

    if plot_mode:
        plt.plot(expanded_time_array, I_1_expanded, label = "time delayed pulse")
        plt.plot(expanded_time_array,I_2_expanded, label = "without random delay")
        plt.legend()
        plt.show()

    

    return I_1_expanded, I_2_expanded, expanded_time_array



# overlap = overlap_percentage(I_1,I_2)



# print("overlap perc: "+str(overlap))



I_1_shifted, I_2_shifted, expanded_time_array = generate_random_phase_delay(I_1,I_2,time_step_fs,time_array_BBO_pulse,plot_mode=True)

# overlap = overlap_percentage(I_1_shifted, I_2_shifted)

# print("overlap perc: "+str(overlap))


def estimate_photon_count(I_1,I_2, time_array):
    mu_0 = 4*np.pi*1e-7  #double check this?
    c = 3e8
    w_plus = 2*np.pi* (1/800 + 1/400)*3e8
    N = 6e19 #need to change this at points!!  #double check this?? N originally in m^-3 then we change to cm^-2?
    beta_plus = 10e-50
    f_rep_rate = 1e3  #in Hz
    h_bar = 6.63e-34 
    length =  5e-6    #size of liquid jet   #NEED TO DOUBLE CHECK THESE!! could shine at an angle to get more length, since would scale with length
    waist = 700e-6   #need to balance this out with intensity!!!  need to use what they give in lightwave explorer at the moment - this will change with using a telescope system! also intensity will change!! - and this system should depend on the spatial thing??
    V = waist*length
    time_integral = np.trapezoid(I_1*I_2, x=time_array)






######now simulate the overlap provided a gradation:


def gradation_check(dx,dt_fs,time_array, I_1, I_2,plot_mode = False):   #assume I_1 is the one we 'move', and I_1 I_2 are after the random phase addition

    dt = time_gradation(dx)*1e15
    # print(dt)
    no_iter = int(np.round((time_array[-1] - time_array[0]) / dt)) + 1
    # print(no_iter)
    time_shift_array = []
    overlap_array = []
    for i in range(-no_iter,no_iter):
        # print(i)
        time_shift = time_gradation(dx*i)*1e15
        time_shift_array.append(time_shift)
        time_delay_array = time_array + time_shift


        t_min = min(time_array[0], time_delay_array[0])
        t_max = max(time_array[-1], time_delay_array[-1])
        num_steps = int(np.round((t_max - t_min) / dt_fs)) + 1
        expanded_time_array = np.linspace(t_min, t_max, num_steps)
        I_1_expanded = np.interp(expanded_time_array, time_delay_array, I_1, left=0.0, right=0.0)   
        I_2_expanded = np.interp(expanded_time_array, time_array, I_2, left=0.0, right=0.0)

        overlap = overlap_percentage(I_1_expanded, I_2_expanded)
        overlap_array.append(overlap)

        # if overlap > np.max(overlap_array):



    if plot_mode:
        print("max achievable efficiency: "+str(max(overlap_array)))
        plt.plot(time_shift_array,overlap_array)
        plt.show()

    return max(overlap_array)




# dx = 5e-6  #10 microns

# gradation_check(dx,time_step_fs,expanded_time_array,I_1_shifted,I_2_shifted)



def efficiency_check(N, I_1, I_2, dt_fs, original_time_array,dx):   #repeat this over several random attempts

    max_array = []
    saved_max = 0
    intensity_product = 0
    # dx = 10e-6
    for i in range(N):

        I_1_expanded, I_2_expanded, expanded_time_array = generate_random_phase_delay(I_1, I_2, dt_fs, original_time_array)

        max_val = gradation_check(dx,dt_fs,expanded_time_array,I_1_expanded,I_2_expanded)

        max_array.append(max_val)

        if max_val > np.max(max_array):

            saved_max = max_val



    
    # plt.plot(max_val)

    mean_val = np.mean(max_array)

    print("mean max efficiency achievable: "+str(mean_val))

    return mean_val, saved_max
   

N = 500
# efficiency_check(N,I_1,I_2,time_step_fs,time_array_BBO_pulse,10e-6)


#do a plot of dx against average max efficiency




def dx_versus_efficiency_plot(dx_array,I_1,I_2):

    mean_val_array = []

    for dx in dx_array: 
        print(dx)

        mean_val = efficiency_check(100,I_1,I_2,time_step_fs,time_array_BBO_pulse,dx)

        mean_val_array.append(mean_val)

    plt.scatter(dx_array*1e6, mean_val_array)   #x,y
    plt.xlabel("dx gradation in microns")
    plt.ylabel("max intensity overlap possible on average")
    plt.show()

dx_array = np.linspace(1,10,100)*1e-6

# dx_versus_efficiency_plot(dx_array, I_1,I_2)


















"""

def align_and_expand_arrays(original_time_array, I_1, I_2, delay_fs):


    # 1. Determine your exact time step (dt)
    dt = original_time_array[1] - original_time_array[0]
    
    # 2. Calculate the shifted coordinates for I_1
    shifted_time_array = original_time_array + delay_fs
    
    # 3. Find the absolute minimum and maximum times needed to hold both pulses
    t_min = min(original_time_array[0], shifted_time_array[0])
    t_max = max(original_time_array[-1], shifted_time_array[-1])
    
    # 4. Create the new Expanded Master Time Array
    # We calculate the exact number of steps to ensure the spacing remains identical
    num_steps = int(np.round((t_max - t_min) / dt)) + 1
    expanded_time_array = np.linspace(t_min, t_max, num_steps)
    
    # 5. Map BOTH pulses onto this new master grid
    # Anything outside a pulse's specific active window is safely zero-padded
    I_1_expanded = np.interp(expanded_time_array, shifted_time_array, I_1, left=0.0, right=0.0)
    I_2_expanded = np.interp(expanded_time_array, original_time_array, I_2, left=0.0, right=0.0)
    
    return expanded_time_array, I_1_expanded, I_2_expanded



"""


















