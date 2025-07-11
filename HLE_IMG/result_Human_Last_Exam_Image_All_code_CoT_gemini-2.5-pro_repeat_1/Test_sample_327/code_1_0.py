# This script classifies three graphene nanoribbon band structures
# based on their visual features and concatenates the results.

# Classification for Ribbon 1
edge_1 = 'A'  # Armchair: parabolic bands with a gap at k=0
width_1 = 8   # Counted 8 bands above the Fermi level (E=0)
band_1 = 1    # Semiconducting: a clear band gap is visible
class_1 = f"{edge_1}{width_1}{band_1}"

# Classification for Ribbon 2
edge_2 = 'A'  # Armchair: bands cross linearly at k=0
width_2 = 7   # Counted 7 bands above the Fermi level (E=0)
band_2 = 0    # Metallic: bands touch at the Fermi level (E=0)
class_2 = f"{edge_2}{width_2}{band_2}"

# Classification for Ribbon 3
edge_3 = 'Z'  # Zigzag: characteristic flat band at E=0
width_3 = 5   # Counted 5 bands above the Fermi level at k=0
band_3 = 0    # Metallic: the flat band at E=0 ensures it's metallic
class_3 = f"{edge_3}{width_3}{band_3}"

# Concatenate the individual classifications to form the final string
final_answer = class_1 + class_2 + class_3

# Print the final concatenated string
print(final_answer)