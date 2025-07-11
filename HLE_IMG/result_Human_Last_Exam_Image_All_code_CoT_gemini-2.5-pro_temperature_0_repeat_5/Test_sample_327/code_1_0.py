# Based on the analysis of the three band structure plots for graphene nanoribbons.

# Classification for Plot 1
edge_1 = 'Z'  # Zigzag edge due to flat bands at Fermi level
width_1 = 8   # 8 bands counted in the conduction band (E > 0)
band_1 = 0    # Metallic as there is no band gap

# Classification for Plot 2
edge_2 = 'A'  # Armchair edge, no flat bands
width_2 = 8   # 8 bands counted in the conduction band
band_2 = 0    # Metallic, bands touch at E=0

# Classification for Plot 3
edge_3 = 'A'  # Armchair edge, no flat bands
width_3 = 5   # 5 bands counted in the conduction band
band_3 = 1    # Semiconducting, clear band gap

# Concatenate the classifications without separators
classification_1 = f"{edge_1}{width_1}{band_1}"
classification_2 = f"{edge_2}{width_2}{band_2}"
classification_3 = f"{edge_3}{width_3}{band_3}"

final_string = classification_1 + classification_2 + classification_3

# The final equation is the concatenation of the individual classifications.
# We print the final result.
print(f"{classification_1}{classification_2}{classification_3}")