# The Chern number is a topological invariant. When two topological systems,
# such as Chern insulators, are brought together without closing the bulk gap,
# the Chern number of the combined system is the sum of the individual Chern numbers.

# Define the Chern number for the first insulator.
c1 = 1

# Define the Chern number for the second insulator.
c2 = 1

# Calculate the Chern number of the junction.
c_total = c1 + c2

# Print the explanation and the final calculation.
print("The Chern number of the junction is the sum of the individual Chern numbers.")
print(f"C_junction = C1 + C2")
print(f"C_junction = {c1} + {c2} = {c_total}")
