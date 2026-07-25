import h5py
import numpy as np
import matplotlib.pyplot as plt

file_path = "my_results/acherus_snapshots.h5"

with h5py.File(file_path,"r") as f:

    # print("Main folders/variables in this file:", list(f.keys()))
    
    # 2. Extract using the exact keys from your terminal output
    e_field_envelope = np.array(f['envelope_snapshot_rzt'])
    plasma_density = np.array(f['density_snapshot_rzt'])
    z_index = np.array(f['snap_z_idx'])

# 3. Check the dimensions!
print(f"Envelope Shape: {e_field_envelope.shape}")
print(f"Density Shape:  {plasma_density.shape}")
print(f"Saved Z-Index:  {z_index}")

# 4. Prepare the 2D slice for plotting
# If the shape is 3D (r, z, t), we want to plot the final z-step.
# We slice it using [:, -1, :] which means (all radial, last z, all time)
if e_field_envelope.ndim == 3:
    plot_data = e_field_envelope[:, 1, :]
    print("Detected a 3D array. Plotting the final z-slice.")
else:
    plot_data = e_field_envelope
    print("Detected a 2D array. Plotting directly.")

# 5. Visualize it
plt.figure(figsize=(8, 6))
# Remember to plot the absolute amplitude of the complex envelope
plt.imshow(np.abs(plot_data), aspect='auto', origin='lower', cmap='inferno')
plt.colorbar(label="E-Field Envelope Amplitude")
plt.title(f"Beam Profile at Snapshot Z-Index: {z_index}")
plt.xlabel("Time (Nodes)")
plt.ylabel("Radial Space (Nodes)")
plt.show()
















