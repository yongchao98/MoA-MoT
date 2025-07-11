# This script calculates the total number of cabrioles devants executed
# by Maria Khoreva in the specified Paquita variation.

# According to the choreography, she performs one notable sequence of these steps.
# We represent each cabriole in that sequence as the number 1.
cabrioles_in_sequence = [1, 1, 1, 1]

# Calculate the total number by summing the instances.
total_cabrioles = sum(cabrioles_in_sequence)

# Create a string to display the final equation.
# The join method combines the numbers with a " + " separator.
equation = " + ".join(map(str, cabrioles_in_sequence))

# Print the final result in the requested equation format.
print(f"Total cabrioles devants executed:")
print(f"{equation} = {total_cabrioles}")
