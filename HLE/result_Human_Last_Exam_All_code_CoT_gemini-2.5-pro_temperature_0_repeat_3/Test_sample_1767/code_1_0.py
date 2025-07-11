# Define the Chern numbers for the two insulators.
c1 = 1
c2 = 1

# The Chern number is an additive topological quantum number. When two
# Chern insulators are joined together with a negligible tunneling barrier,
# the resulting system can be described by a total Chern number that is
# the sum of the individual Chern numbers.
# The counter-propagating edge states at the interface gap out,
# leaving a new bulk system with a combined topological charge.

# Calculate the total Chern number of the junction.
c_total = c1 + c2

# Print the explanation and the final calculation.
print("The Chern number of the combined system is the sum of the individual Chern numbers.")
print(f"Chern number of the first insulator (C1): {c1}")
print(f"Chern number of the second insulator (C2): {c2}")
print("\nThe total Chern number (C_total) is C1 + C2.")
print("The final equation is:")
print(f"{c_total} = {c1} + {c2}")