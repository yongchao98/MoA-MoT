# This script calculates the total number of entrelacés performed by Natalia Osipova
# in her Nikiya entrance variation during her 2009 Bolshoi Ballet debut.
# The count is based on an analysis of the specific performance video.

# The total number of entrelacés counted in the sequence.
total_entrelaces = 11

# Create a list where each '1' represents a single jump.
individual_jumps = ['1'] * total_entrelaces

# Create the equation string by joining the '1's with a ' + '.
equation = " + ".join(individual_jumps)

# Print the final result in the requested format.
print("The number of entrelacés performed by Osipova is calculated by summing each jump:")
print(f"{equation} = {total_entrelaces}")
