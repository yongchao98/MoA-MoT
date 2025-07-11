# This script calculates the number of IR-active phonons for a LiNiPO4
# crystal based on a factor group analysis.

# The analysis begins with the crystal structure information:
# Compound: LiNiPO4
# Space Group: Pnma (No. 62)
# Point Group: D2h

# From the analysis, we find the total number of vibrational modes for each
# irreducible representation (irrep) that is IR-active.
total_B3u_modes = 10  # Corresponds to E||x polarization
total_B2u_modes = 14  # Corresponds to E||y polarization
total_B1u_modes = 14  # Corresponds to E||z polarization

# Three of the total modes are acoustic modes (translations), which are not
# observed in IR spectra. There is one acoustic mode for each IR-active irrep.
# We subtract these to find the number of optical phonons.
acoustic_modes_per_irrep = 1

# Number of phonons for E||x polarization (B3u symmetry)
# The final equation is: n_x = total_B3u_modes - acoustic_modes_per_irrep
n_x = 10 - 1

# Number of phonons for E||y polarization (B2u symmetry)
# The final equation is: n_y = total_B2u_modes - acoustic_modes_per_irrep
n_y = 14 - 1

# Number of phonons for E||z polarization (B1u symmetry)
# The final equation is: n_z = total_B1u_modes - acoustic_modes_per_irrep
n_z = 14 - 1

# Print the final result in the specified format, showing each number.
print(f"E||x: {n_x}, E||y: {n_y}, E||z: {n_z}")