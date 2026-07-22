
from __future__ import annotations
import os
from pathlib import Path
import numpy as np

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

    print("SIZE of laser field: ")
    print(laser_field.shape)



    # =====================================================================
    # 4. THE INJECTION POINT
    # =====================================================================


    # print("Loading custom wavefront...")
    # my_custom_array = np.load("pre_curved_wavefront.npy")
    
    ## Overwrite the built-in Gaussian array with your generic beam!
    ## (Note: Check if the attribute is .field, .E, or .envelope in laser.py)
    # laser.init_envelope = lambda: my_custom_array
 




 

    # # 5. RESUME NORMAL PHYSICS: Now that your array is in, build the equations
    # equation = Equation(medium, laser, grid)
    # ionization = KeldyshIonization(
    #     medium, laser, model_name=config.ionization_model, params=config.ionization_model_par
    # )

    # # 6. RUN THE SOLVER
    # print("Starting propagation...")
    # solver = init_solver(config, medium, laser, grid, equation, ionization, output)
    # fft_manager.set_fft_backend(config.computing_backend)
    
    # solver.propagate()
    # output.save_results(solver, grid)
    # print("Simulation complete! Results saved.")

if __name__ == "__main__":
    run_my_simulation()












































