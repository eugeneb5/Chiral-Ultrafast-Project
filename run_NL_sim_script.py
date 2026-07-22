
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
# import tomllib
# from pathlib import Path

# Import the main setup functions and classes from Acherus
from acherus.__main__ import load_config, init_solver
from acherus.mesh.grid import Grid
from acherus.physics.medium import Medium
from acherus.physics.laser import Laser
from acherus.physics.equation import Equation
from acherus.physics.keldysh import KeldyshIonization
from acherus.data.store import OutputManager
from acherus.functions.fft_backend import fft_manager



def run_my_simulation():
    # 1. THE SETUP: Let Acherus read your TOML file to get the grid sizes
    config_path = Path("002_air_IR_fs.toml")
    config = load_config(config_path)
    
    output_dir = Path("./my_results")
    os.environ["ACHERUS_BASE_DIR"] = str(output_dir)    #DOUBLE CHECK THIS BIT??
    output = OutputManager(save_path=output_dir)




    # 2. BUILD THE UNIVERSE: Create the empty space and air medium
    grid = Grid(space_par=config.space_par, axis_par=config.axis_par, time_par=config.time_par)
    medium = Medium(medium_name=config.medium_name, medium_par=config.medium_par)
    





    # 3. BUILD THE LASER: This creates the default Gaussian pulse in memory
    laser = Laser(medium, grid, pulse_name=config.pulse_name, pulse_par=config.pulse_par)


    laser_field = laser.init_envelope()

    # print("SIZE of laser field: ")
    # print(laser_field.shape)

    # plt.imshow(np.abs(laser_field)**2, origin='lower', aspect='auto', cmap='Blues')
    # plt.xlabel("time")
    # plt.ylabel("x")
    # plt.show()


    # time_centre = laser_field.shape[1]//2
    # plt.figure(figsize=(10, 6))
    # plt.imshow(np.abs(laser_field)**2, aspect='auto', origin='lower', cmap='Blues') 
    # plt.axvline(x=time_centre, color='red', linestyle='--', linewidth=2, 
    #                 label=f'Proposed Chop Line (t = {time_centre})')
    # plt.show()


    # =====================================================================
    # 4. THE INJECTION POINT
    # =====================================================================


    print("Loading custom wavefront...")
    my_custom_array = np.load("E_field_z_0.npy")
    fixed_array = my_custom_array.T.copy()

 

    # center_index = fixed_array.shape[0] // 2
    # amplitude = np.abs(fixed_array)
    # plt.figure(figsize=(10, 6))
    # plt.imshow(amplitude, aspect='auto', origin='lower', cmap='Blues') 
    # plt.axhline(y=center_index, color='red', linestyle='--', linewidth=2, 
    #             label=f'Proposed Chop Line (x = {center_index})')

    # plt.colorbar(label="E-Field Amplitude")
    # plt.title("Verification Check: Does the mathematical middle match the beam peak?")
    # plt.xlabel("Time (Nodes)")
    # plt.ylabel("Spatial Axis (Nodes)")
    # plt.legend()
    # plt.show()
 


    center_index = fixed_array.shape[0] // 2
    radial_array = fixed_array[center_index:, :]
    laser.init_envelope = lambda: radial_array  #replace function in class!!!


    laser_field = laser.init_envelope()
    # plt.imshow(np.abs(laser_field)**2, origin='lower', aspect='auto', cmap='Blues')
    # plt.xlabel("time")
    # plt.ylabel("x")
    # plt.show()



    #check if goes thru centre of time part?!


    time_centre = laser_field.shape[1]//2

    # plt.figure(figsize=(10, 6))
    # plt.imshow(np.abs(laser_field)**2, aspect='auto', origin='lower', cmap='Blues') 
    # plt.axvline(x=time_centre, color='red', linestyle='--', linewidth=2, 
    #                 label=f'Proposed Chop Line (t = {time_centre})')
    # plt.show()


    on_axis_field = np.abs(laser_field[0,:])
    # plt.plot(on_axis_field)
    # plt.show()
    current_peak = np.argmax(on_axis_field)
    shift_amount = time_centre - current_peak
    centred_laser_field = np.roll(laser_field, shift_amount, axis=1)





    laser.init_envelope = lambda: centred_laser_field  #replace function in class!!!


    check_var = laser.init_envelope()

    plt.figure(figsize=(10, 6))
    plt.imshow(np.abs(check_var)**2, aspect='auto', origin='lower', cmap='Blues') 
    plt.axvline(x=time_centre, color='red', linestyle='--', linewidth=2, 
                    label=f'midpoint line (t = {time_centre})')


    plt.show()


















    # 5. RESUME NORMAL PHYSICS: Now that your array is in, build the equations
    equation = Equation(medium, laser, grid)
    ionization = KeldyshIonization(
        medium, laser, model_name=config.ionization_model, params=config.ionization_model_par
    )

    # 6. RUN THE SOLVER
    print("Starting propagation...")
    solver = init_solver(config, medium, laser, grid, equation, ionization, output)
    fft_manager.set_fft_backend(config.computing_backend)
    
    solver.propagate()
    output.save_results(solver, grid)
    print("Simulation complete! Results saved.")

if __name__ == "__main__":
    run_my_simulation()












































