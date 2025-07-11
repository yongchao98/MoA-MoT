# The total number of modes for each IR-active symmetry in LiNiPO4 (Pnma space group)
# These values are derived from a full factor group analysis.
# B1u corresponds to E||z, B2u to E||y, and B3u to E||x.
total_modes_B3u = 14 # for E||x polarization
total_modes_B2u = 14 # for E||y polarization
total_modes_B1u = 10 # for E||z polarization

# Each IR-active symmetry has one acoustic mode, which is not observed in the optical spectrum.
acoustic_modes_per_symmetry = 1

# Calculate the number of observable optical phonons for each polarization
# by subtracting the acoustic mode.
phonons_x = total_modes_B3u - acoustic_modes_per_symmetry
phonons_y = total_modes_B2u - acoustic_modes_per_symmetry
phonons_z = total_modes_B1u - acoustic_modes_per_symmetry

# Print the final result in the requested format
print(f"E||x: {phonons_x}, E||y: {phonons_y}, E||z: {phonons_z}")
