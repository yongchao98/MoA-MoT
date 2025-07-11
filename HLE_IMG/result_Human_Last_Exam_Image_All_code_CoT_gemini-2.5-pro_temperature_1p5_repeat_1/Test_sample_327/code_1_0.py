# Step 1: Define the classifications for each nanoribbon based on visual analysis.

# Analysis for Plot 1
# Edge Type: Armchair (A), as band gap is at k=0.
# Width (N): 8 bands are visible in the conduction band (E > 0).
# Band Type: Metallic (0), as the valence and conduction bands touch at E=0.
edge1 = "A"
width1 = 8
band1 = 0

# Analysis for Plot 2
# Edge Type: Zigzag (Z), due to the characteristic flat edge state bands near E=0.
# Width (N): 7 bands are visible in the conduction band at k=0.
# Band Type: Metallic (0), as is always the case for ideal zigzag nanoribbons.
edge2 = "Z"
width2 = 7
band2 = 0

# Analysis for Plot 3
# Edge Type: Armchair (A), as band gap is at k=0.
# Width (N): 7 bands are visible in the conduction band.
# Band Type: Semiconducting (1), as there is a clear band gap at E=0.
edge3 = "A"
width3 = 7
band3 = 1

# Step 2: Concatenate the classifications into the final required format "EdgeWidthBand".
part1 = f"{edge1}{width1}{band1}"
part2 = f"{edge2}{width2}{band2}"
part3 = f"{edge3}{width3}{band3}"

final_string = f"{part1}{part2}{part3}"

# Step 3: Print the final concatenated string.
# The following line can be interpreted as the "final equation" mentioned in the instructions,
# showing how the final string is constructed from its parts.
print(f"Classification 1: {edge1}{width1}{band1}")
print(f"Classification 2: {edge2}{width2}{band2}")
print(f"Classification 3: {edge3}{width3}{band3}")
print(f"Final Concatenated String: {final_string}")
