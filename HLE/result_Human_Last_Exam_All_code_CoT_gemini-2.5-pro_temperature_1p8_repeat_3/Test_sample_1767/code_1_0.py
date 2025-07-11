# Define the Chern numbers of the two insulators
c1 = 1
c2 = 1

# The Chern number of a junction is the sum of the individual Chern numbers
c_junction = c1 + c2

# Print the explanation and the final equation
print("The Chern number of a system formed by joining two topological materials is the sum of their individual Chern numbers.")
print(f"Chern number of the first insulator: C1 = {c1}")
print(f"Chern number of the second insulator: C2 = {c2}")
print("The Chern number of the junction is C_junction = C1 + C2.")
print(f"C_junction = {c1} + {c2} = {c_junction}")