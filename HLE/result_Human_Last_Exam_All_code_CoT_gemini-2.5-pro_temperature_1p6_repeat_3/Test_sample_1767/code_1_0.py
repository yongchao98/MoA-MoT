# Define the Chern numbers for the two insulators
c1 = 1
c2 = 1

# The Chern number of the junction is the sum of the individual Chern numbers
# because the Chern number is an additive topological invariant.
c_total = c1 + c2

# Print the explanation and the final calculation
print(f"The first Chern insulator has a Chern number C1 = {c1}.")
print(f"The second Chern insulator has a Chern number C2 = {c2}.")
print("The Chern number of the combined junction is the sum of the individual numbers.")
print(f"C_total = C1 + C2")
print(f"C_total = {c1} + {c2} = {c_total}")