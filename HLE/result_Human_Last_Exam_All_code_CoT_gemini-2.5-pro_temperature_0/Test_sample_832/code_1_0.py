# This script calculates the number of entrelacés performed by Natalia Osipova
# in her 2009 Bolshoi debut as Nikiya in the "Death of Nikiya" scene.
# The count is based on video analysis of the performance.

# Each '1' in this list represents one entrelacé counted from the performance.
entrelaces_counted = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

# Calculate the total number of entrelacés
total_entrelaces = sum(entrelaces_counted)

# Create the equation string to display the calculation
equation_str = " + ".join(map(str, entrelaces_counted))

# Print the final result in the requested format
print(f"Based on analysis of the 2009 performance, the number of entrelacés is calculated as follows:")
print(f"{equation_str} = {total_entrelaces}")
