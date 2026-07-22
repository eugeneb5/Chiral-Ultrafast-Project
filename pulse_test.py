import numpy as np
import matplotlib.pyplot as plt
import LightwaveExplorer as Le
from scipy.signal import hilbert

# help(Le)

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
    time_span = sim.timeSpan*1e15 #in fs
    time_step = sim.timeStep  #is in s
    time_step_fs = time_step 
    
    print(N_time, time_span,time_step)

    spatial_step = sim.spatialStep
    # spatial_height = sim.spatialHeight  #this is y var - so ignore since used 2D simulation!
    spatial_width = sim.spatialWidth* 1e6  #this is x var - so use this! in microns here
    N_space = sim.Nspace

    print(N_space,spatial_width )

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


    return N_time, time_step_fs, full_E_field, space_array*1e-6




results_zip = "Data_BBO_sim/BBO_USE_THIS_TEST.zip"

no_BBO_pulse = "Data_BBO_sim/New folder/test_1.zip"

test_pulse_2 = "Data_BBO_sim/New folder/test_1_400nm.zip"

acherus_test = "test_pulse_for_acherus_1.zip"





# N_time,  time_step_fs ,full_E_field_BBO ,space_array = define_pulse(results_zip, read_char = False)   #space_array is in microns!

N_time, time_step_fs,full_E_field_no_BBO,space_array = define_pulse(acherus_test,read_char = False ,pol_y=True)

# N_time_2, time_step_fs_2,E_field_2,space_array_2 = define_pulse(test_pulse_2,read_char = False ,pol_y=False)








# E_xt_no_BBO = hilbert(full_E_field_no_BBO, axis=0)

# plt.imshow(np.abs(E_xt_no_BBO)**2, origin='lower', aspect='auto')
# plt.show()



plt.imshow(np.abs(hilbert(full_E_field_no_BBO, axis=0))**2, origin='lower', aspect='auto', cmap='Blues')

# plt.imshow(np.abs(E_xt_no_BBO)**2, origin='lower', aspect='auto', cmap='Reds', alpha=0.5)

# plt.imshow(np.abs(hilbert(E_field_2, axis=0))**2, origin='lower', aspect='auto', cmap='Reds',alpha = 0.5)

plt.xlabel("Space (x) / microns")
plt.ylabel("Time (t) / fs")

plt.show()
















def time_to_angular_freq_space(U_xt,time_step_fs,N_time,space_array):

    # print(U_xt.shape)
    U_xt_complex = hilbert(U_xt, axis=0) #to get rid of negative frequency parts
    N_t = N_time
    dt = time_step_fs
    dw = 2*np.pi/(N_t*dt)
    w_max = 2*np.pi/dt  #so -w_max/2 to w_max/2 !!! this is just the full range, so 0 to w_max/2 are the positive frequencies, therefore need smaller spacings for more w values!! makes sense...

    plot_extent = [space_array[0],space_array[-1],-w_max/2,w_max/2,]   #space array is in microns!

    # w_array = np.linspace(-w_max/2, w_max/2 +1, N_t)

    w_array = np.fft.fftshift(np.fft.fftfreq(U_xt.shape[0], d=dt)) * 2 * np.pi



    # w_min = w_array[0]
    # w_max = w_array[-1]
    # plot_extent = [space_array[0], space_array[-1], w_min, w_max]





    # plt.imshow(U_xt, origin='lower', aspect='auto', cmap='viridis')
    # plt.show()
    U_xw_raw = np.fft.fft(U_xt_complex, axis=0)   #axis 0 is time axis...
    U_xw = np.fft.fftshift(U_xw_raw, axes=0)
    plt.imshow(np.abs(U_xw), extent =plot_extent,origin='lower', aspect='auto', cmap='viridis')
    # plt.ylim(2.2e15, 2.5e15)  #range much smaller than where its centred(which is based on carrier freq!!)
    plt.xlabel("x /microns")
    plt.ylabel("w / rads^-1")
    plt.show()
    
    # print(U_xw.shape)

    return U_xw, w_array





U_xw, w_array = time_to_angular_freq_space(full_E_field_no_BBO,time_step_fs,N_time,space_array)

# U_xw_2, w_array_2 = time_to_angular_freq_space(E_field_2, time_step_fs_2,N_time_2, space_array_2 )







def ref_index(url, w):  #for single w value...
    sellmeier_coeffs = Le.getSellmeierFromRII(url)

    wavelength = 2*np.pi*3e8 / w *1e6   #need it in microns!

    wavelength_array = np.array([wavelength])

    n = Le.sellmeier(wavelength_array, sellmeier_coeffs,0) 

    return n[0]


def ref_index_array(url, w_array):  
    sellmeier_coeffs = Le.getSellmeierFromRII(url)

    wavelength_array = 2*np.pi*3e8 / w_array *1e6

    n_array = Le.sellmeier(wavelength_array, sellmeier_coeffs,0) 

    return n_array


def ASW_propogation(U_x_slice, z, space_array, w_i):   #archaic function, dont need



    N_x = len(space_array)
    dx = space_array[1] - space_array[0]
    kx_array = np.fft.fftshift(np.fft.fftfreq(N_x, d=dx)) * 2 * np.pi
    k = w_i/3e8  #k_magnitude!


    U_kx = np.fft.fftshift(np.fft.fft(U_x_slice))

    propogation_factor = np.exp(1j * k * z) * np.exp(-1j * (kx_array**2) / (2 * k) * z)
    U_kx_propagated = U_kx * propogation_factor

    U_slice_z = np.fft.ifft(np.fft.ifftshift(U_kx_propagated))

    return U_slice_z


def Lens_diameter(R,L_0):

    return 2*np.sqrt(R**2 - (R-L_0)**2)



def lens_and_propagation(U_xw, z, material_url,R_lens,D_lens, space_array=space_array, w_array=w_array):  #space_array is in microns!!! but please have R and D in terms of metres!!!!


    U_xw_after_lens = np.zeros_like(U_xw, dtype=complex)
    c = 3e8
    # space_array_m = space_array #convert back into metres here!

    def L_x(R,D,x):

        sum_1 = np.sqrt(R**2 - x**2)
        sum_2 = np.sqrt(R**2 - (D/2)**2)
        
        sum = sum_1 - sum_2

        # sum = L_0 - (R - np.sqrt(R**2 - x**2))

        assert np.all(np.isreal(sum)), "some elements are complex!"

        return sum
    
    n_array = ref_index_array(material_url, w_array) #should be same size as w_array... could check?

    assert len(n_array) == len(w_array), "something wrong with ref index array" 

    N_x = len(space_array)
    dx = space_array[1] - space_array[0]
    # print("checking values")
    # print(dx, N_x)
    # print(space_array)
    kx_array = np.fft.fftshift(np.fft.fftfreq(N_x, d=dx)) * 2 * np.pi  #is the same for all of w_i, x arrays - so just do once...


    L = L_x(R_lens, D_lens, space_array)  #only need once...

    U_w_size = U_xw.shape[0]



    for i in range(U_w_size):  #this elects w dim right?
        
        if w_array[i]<=0:
            continue  #i.e. restart if is a negative freq..
       
        # print("iteration index: "+str(i) +"of max index" +str(U_w_size))
        w_i = w_array[i]
        k = w_i/3e8  #k_magnitude!
        # wavelength = 2*np.pi/k *1e9  #so in nm!
        n = n_array[i]
        # print(w_i,n,wavelength)
        U_xwi = U_xw[i,:]  #this is the w_i slice
        

        #phase mask first
        phi = (w_i/c)*(n-1)*L

        #apply lens phase effect

        U_xwi_after_lens = U_xwi * np.exp(1j * phi)

        #now propogate using ASW
        # U_kx = np.fft.fftshift(np.fft.fft(U_xwi_after_lens))
        U_kx = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(U_xwi_after_lens)))
        # propogation_factor = np.exp(1j * k * z) * np.exp(-1j * (kx_array**2) / (2 * k) * z)
        propogation_factor = np.exp(-1j * (kx_array**2) / (2 * k) * z)
        U_kx_propagated = U_kx * propogation_factor

        #save array
        # U_xw_after_lens[i, :] = np.fft.ifft(np.fft.ifftshift(U_kx_propagated))
        U_xw_after_lens[i, :] = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(U_kx_propagated)))


    #now transform back to time domain...
    U_xw_unaligned = np.fft.ifftshift(U_xw_after_lens, axes=0)
    U_xt_final_complex = np.fft.ifft(U_xw_unaligned, axis=0)
    E_field_final_real = 2*np.real(U_xt_final_complex)  #double check if need???


    plt.imshow(np.abs(U_xt_final_complex)**2, origin='lower', aspect='auto', cmap='viridis')
    plt.xlabel("Space (x) / microns")
    plt.ylabel("Time (t) / fs")
    plt.show()

    return E_field_final_real



#5mm waist at lens for this test...

def Z_R_estimate(wavelength, waist_lens, focal_length):

    return (focal_length**2) * wavelength / (np.pi*waist_lens**2)


###want to make something that will tell us if Z_R >> liquid jet thickness, then can shortcut the process for calculating intensity...


def save_array(save_name,E_array):   #name; #file saved as "..., .npy"



    np.save(save_name, E_array)
    print("Data saved")


def load_array(load_name, plot_mode = False):   #file saved as "..., .npy"



    U_xt = np.load(load_name)

    

    if plot_mode:

        plt.imshow(np.abs(hilbert(U_xt,axis=0))**2, origin='lower', aspect='auto')
        plt.xlabel("Space (x) / microns")
        plt.ylabel("Time (t) / fs")
        plt.show()

    return U_xt
    




def temporal_shift(U_xt, U_xt_2, time_array):   #work in progress...


    I_1 = np.abs(hilbert(U_xt, axis=0))**2
    I_2 = np.abs(hilbert(U_xt_2, axis = 0))**2

    total_time_pixels = I_1.shape[0]

   
    max_shift = total_time_pixels // 2

 
    pixel_shifts = np.arange(-max_shift, max_shift + 1, 1)



    r_coords = np.abs(space_array) * 1e-6  # Convert to meters
    _, R_matrix = np.meshgrid(time_array, r_coords, indexing='ij')
    volume_weighting = 2 * np.pi * R_matrix

    signal_curve = []

    

def focal_length_at_wavelength(R,url,wavelength):   #for a planoconvex lens!
    w = 3e8*2*np.pi/wavelength
    n_ = np.real(ref_index(url, w))
    print("ref index retrieved: "+str(n_))
    
    return R/(n_-1)



def round_to_nearest_mm(value_in_meters):
    """
    Rounds a meter value to the nearest millimeter (3 decimal places).
    The f-string formatting {:.3f} ensures that numbers ending in zero 
    (like 0.18) are displayed as 0.180 instead of being truncated.
    """

    rounded_value = round(value_in_meters, 3)
    return rounded_value




def round_to_nearest_cm(value_in_meters):
    """
    Rounds a meter value to the nearest centimeter (2 decimal places),
    but pads the string to 3 decimal places to display the mm placeholder.
    This turns 0.179 into 0.180.
    """
    rounded_value = round(value_in_meters, 2)
    return rounded_value









# run simulation!

# jet_width = 50e-6
# z_array = [focal_point- jet_width/2, focal_point, focal_point+jet_width/2]

url = "https://refractiveindex.info/database/data/main/MgF2/nk/Li-o.yml"
R_lens = 75.4e-3  #0.36
L_0 = 3.2e-3  #93e-3
D_lens = 25.4e-3

focal_point = focal_length_at_wavelength(R_lens, url, 810e-9)  #expect in mm

z_0 = round_to_nearest_mm(focal_point)

print(focal_point, z_0)



z_array = [0,z_0/2 , z_0, 3*z_0/2]


for z in z_array:

    print("starting calculation for z = "+str(z))

    E_field = lens_and_propagation(U_xw, z, url, R_lens,D_lens, space_array = space_array, w_array = w_array)
    
    # E_field = lens_and_propagation(U_xw_2, z, url, R_lens,D_lens, space_array = space_array_2, w_array = w_array_2)


    save_array(f"E_field_z_{z}.npy", E_field)





U_xt_after = load_array("E_field_z_0.npy",plot_mode= True)

U_xt_after_2 = load_array("E_field_z_0.1005.npy",plot_mode= True)



# plt.imshow(np.abs(hilbert(U_xt_after, axis=0))**2, origin='lower', aspect='auto', cmap='Blues')

# # plt.imshow(np.abs(E_xt_no_BBO)**2, origin='lower', aspect='auto', cmap='Reds', alpha=0.5)

# plt.imshow(np.abs(hilbert(U_xt_after_2, axis=0))**2, origin='lower', aspect='auto', cmap='Reds',alpha = 0.5)

# plt.xlabel("Space (x) / microns")
# plt.ylabel("Time (t) / fs")

# plt.show()


































