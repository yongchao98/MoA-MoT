# Number of IR-active phonons for LiNiPO4 based on factor group analysis.
# The space group is Pnma (D_2h^16).
# The standard crystallographic setting is a, b, c.
# We assume the polarization directions x, y, z correspond to a, b, c respectively.

# The IR-active modes are decomposed as: Gamma_IR = 11*B3u + 11*B1u + 8*B2u
# B3u corresponds to polarization along the a-axis (x)
# B1u corresponds to polarization along the b-axis (y)
# B2u corresponds to polarization along the c-axis (z)

num_phonons_x = 11
num_phonons_y = 11
num_phonons_z = 8

# Print the result in the requested format
print(f"E||x: {num_phonons_x}, E||y: {num_phonons_y}, E||z: {num_phonons_z}")