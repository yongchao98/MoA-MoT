# The Chern number is a topological invariant that characterizes the bulk of a 2D material.
# For non-interacting or weakly coupled systems, the Chern number is an additive quantity.
# When two systems are brought together to form a junction, the Chern number of the
# combined system is the sum of the individual Chern numbers.

# Define the Chern numbers of the two insulators.
c1 = 1
c2 = 1

# Calculate the Chern number of the junction.
c_total = c1 + c2

# Print the explanation and the result.
print("Principle: The Chern number of a combined system is the sum of the Chern numbers of its individual components.")
print(f"Chern number of the first insulator (C1): {c1}")
print(f"Chern number of the second insulator (C2): {c2}")
print("\nThe Chern number of the junction (C_total) is calculated as:")
print(f"C_total = C1 + C2")
print(f"C_total = {c1} + {c2} = {c_total}")
print(f"\nTherefore, the Chern number of the junction is {c_total}.")
