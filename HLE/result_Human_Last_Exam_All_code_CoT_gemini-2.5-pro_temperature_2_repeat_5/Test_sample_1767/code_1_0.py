# Define the Chern numbers for the two insulators
c1 = 1
c2 = 1

# The Chern number of the combined system (the junction) is the sum
# of the individual Chern numbers because the Chern number is an additive topological invariant.
c_junction = c1 + c2

# Print the calculation steps and the final result
print("The Chern number of the junction is the sum of the individual Chern numbers.")
print(f"C_junction = {c1} + {c2}")
print(f"C_junction = {c_junction}")