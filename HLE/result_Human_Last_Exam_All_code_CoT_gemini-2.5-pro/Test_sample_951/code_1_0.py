# Step 1: Define the total number of vibrational modes for each IR-active symmetry.
# These numbers are derived from a group theory analysis of the LiNiPO4 crystal structure (Pnma space group).
total_modes_B3u = 12  # Symmetry corresponding to x-polarization
total_modes_B2u = 9   # Symmetry corresponding to y-polarization
total_modes_B1u = 12  # Symmetry corresponding to z-polarization

# Step 2: The 3 acoustic modes have B1u, B2u, and B3u symmetries.
# Subtract one acoustic mode from each total to find the number of IR-active optical phonons.
num_phonons_x = total_modes_B3u - 1
num_phonons_y = total_modes_B2u - 1
num_phonons_z = total_modes_B1u - 1

# Step 3: Print the final result in the requested format.
# The final equation consists of the number of phonons for each polarization.
print(f"E||x: {num_phonons_x}, E||y: {num_phonons_y}, E||z: {num_phonons_z}")