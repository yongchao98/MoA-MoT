# This script calculates the number of entrelacés in Natalia Osipova's
# "Death of Nikiya" variation from her 2009 Bolshoi debut.
# The number is based on direct observation of the performance recording.

# Each '1' in the list represents one entrelacé performed.
entrelaces_count = [1, 1, 1, 1, 1, 1, 1, 1]

# Calculate the total number of entrelacés
total = sum(entrelaces_count)

# Create the equation string for printing
equation = " + ".join(map(str, entrelaces_count))

# Print the final result in the required format
print("The total number of entrelacés performed is the sum of each individual jump:")
print(f"{equation} = {total}")
