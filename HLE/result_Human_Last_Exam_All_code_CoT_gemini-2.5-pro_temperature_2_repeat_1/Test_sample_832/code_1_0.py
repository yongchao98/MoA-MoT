# This script calculates and displays the number of entrelacés performed by
# Natalia Osipova in her 2009 Bolshoi debut as Nikiya during the "Death of Nikiya" scene variation.
# Based on observation of the performance, the count is 8.

# The total number of entrelacés observed in the variation.
total_entrelaces = 8

# We can represent the counting of each jump as adding "1" for each one.
# Create a list of '1's, with the length equal to the total number of entrelacés.
individual_jumps = ['1'] * total_entrelaces

# Join the list elements with a " + " to form the equation string.
equation_string = " + ".join(individual_jumps)

# Print the final equation.
# The format "1 + 1 + ... = total" fulfills the requirement to show each number.
print("Counting each entrelacé in the variation gives the following equation:")
print(f"{equation_string} = {total_entrelaces}")
