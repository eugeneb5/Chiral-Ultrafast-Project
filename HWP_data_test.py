import numpy as np
import matplotlib.pyplot as plt


filename = 'HWP_dataset_1.txt'

def read_angle_power(filename, plot_mode=False):

    HWP_angle = []
    power = []


    with open(filename, 'r') as file:
    
        next(file) 
        
        for line in file:
            # Strip invisible newline characters and split by the comma
            parts = line.strip().split(',')
            
            # Convert the text into decimal numbers (floats) and add to our arrays
            HWP_angle.append(float(parts[0]))
            power.append(float(parts[1]))
   
    if plot_mode:


        print("Angles:", HWP_angle)
        print("Powers:", power)



        plt.plot(HWP_angle, power)
      
        plt.show()



    return HWP_angle, power


read_angle_power(filename, plot_mode = True)


def read_time_power(filename, plot_mode = False):

    
    power = []

    # Open the text file
    with open(filename, 'r') as file:
        # Skip the "angle" header row
        next(file) 
        
        for line in file:
            # Clean the line and convert directly to a decimal number
            power.append(float(line.strip()))

    array_size = len(power)
    time_array = np.linspace(0,(array_size-1)*5, array_size)

    if plot_mode:
        print("power List:")
        print(power)
        print()
        print("TIME LIST:")
        print(time_array)

        plt.plot(time_array, power)
        plt.show()

    return power, time_array

filename_2 = 'HWP_dataset_2.txt'

# read_time_power(filename_2, plot_mode=True)

   
   
















