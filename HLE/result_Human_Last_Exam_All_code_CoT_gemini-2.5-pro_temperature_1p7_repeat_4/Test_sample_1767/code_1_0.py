# The Chern number is a topological invariant. When two topological systems,
# such as Chern insulators, are joined together with negligible interaction,
# the Chern number of the combined system is the sum of the individual Chern numbers.

# Define the Chern number for each of the two insulators.
c1 = 1
c2 = 1

# Calculate the Chern number of the junction by adding the individual Chern numbers.
c_total = c1 + c2

# Print the formula and the final calculation.
print("The Chern number of the junction is the sum of the individual Chern numbers:")
print(f"C_junction = C1 + C2")
print(f"C_junction = {c1} + {c2}")
print(f"The total Chern number of the junction is: {c_total}")
