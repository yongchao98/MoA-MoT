# Based on the analysis of the three band structure plots.

# Plot 1 classification
edge1 = 'A'  # Armchair, because it is semiconducting
width1 = 8   # Number of conduction bands
band1 = 1    # Semiconducting (has a band gap)
class1 = f"{edge1}{width1}{band1}"

# Plot 2 classification
edge2 = 'Z'  # Zigzag, because it's metallic with N=7 (while N=7 Armchair is semiconducting)
width2 = 7   # Number of conduction bands
band2 = 0    # Metallic (no band gap)
class2 = f"{edge2}{width2}{band2}"

# Plot 3 classification
edge3 = 'A'  # Armchair, because it is semiconducting
width3 = 7   # Number of conduction bands
band3 = 1    # Semiconducting (has a band gap)
class3 = f"{edge3}{width3}{band3}"

# Concatenate the classifications
final_string = class1 + class2 + class3

# Print the breakdown and the final result
# The prompt asks to "output each number in the final equation"
# Here we show the components that form the final string.
print(f"Plot 1: Edge='{edge1}', Width={width1}, Band={band1} -> {class1}")
print(f"Plot 2: Edge='{edge2}', Width={width2}, Band={band2} -> {class2}")
print(f"Plot 3: Edge='{edge3}', Width={width3}, Band={band3} -> {class3}")
print("Final concatenated classification:")
print(final_string)
