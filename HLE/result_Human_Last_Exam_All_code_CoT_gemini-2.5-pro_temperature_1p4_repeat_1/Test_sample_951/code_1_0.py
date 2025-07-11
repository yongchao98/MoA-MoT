# Based on group theory for LiNiPO4 (orthorhombic, Pnma space group).
# Total modes are decomposed into irreducible representations of the D2h point group.
# Gamma_total = 12Ag + 9B1g + 12B2g + 9B3g + 9Au + 12B1u + 9B2u + 12B3u

# Total number of modes for each IR-active symmetry
total_B1u = 12
total_B2u = 9
total_B3u = 12

# There is one acoustic mode for each of these symmetries, which must be subtracted.
acoustic_modes = 1

# Calculate the number of IR-active optical phonons for each symmetry.
# The calculation for each symmetry is explicitly shown below.

# B1u modes correspond to polarization E || c-axis (z-axis)
num_phonons_B1u = total_B1u - acoustic_modes
print(f"Calculation for E||z (B1u modes): {total_B1u} (total) - {acoustic_modes} (acoustic) = {num_phonons_B1u} phonons")

# B2u modes correspond to polarization E || b-axis (y-axis)
num_phonons_B2u = total_B2u - acoustic_modes
print(f"Calculation for E||y (B2u modes): {total_B2u} (total) - {acoustic_modes} (acoustic) = {num_phonons_B2u} phonons")

# B3u modes correspond to polarization E || a-axis (x-axis)
num_phonons_B3u = total_B3u - acoustic_modes
print(f"Calculation for E||x (B3u modes): {total_B3u} (total) - {acoustic_modes} (acoustic) = {num_phonons_B3u} phonons")

# Assigning to x, y, z assuming a standard axis correspondence (x||a, y||b, z||c)
num_x = num_phonons_B3u
num_y = num_phonons_B2u
num_z = num_phonons_B1u

# Print the final result in the requested format
print("\nFinal Answer:")
print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")
