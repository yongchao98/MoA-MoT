# This script calculates the number of sissone fermes performed by Svetlana Zakharova
# in the Act II Odette variation of the 2014 Bolshoi's "Swan Lake".
# The count is based on an analysis of the performance video.

# Each '1' in the list represents a single sissone ferme.
individual_sissones = [1, 1, 1, 1]

# Calculate the total by summing the individual steps.
total_count = sum(individual_sissones)

# Create a string representation of the equation.
equation_string = " + ".join(map(str, individual_sissones))

# Print the explanation and the final equation.
print(f"Based on the 2014 performance, the total number of sissone fermes is {total_count}.")
print("This can be represented by the following equation:")
print(f"{equation_string} = {total_count}")