# Define the Chern numbers of the two individual insulators.
C1 = 1
C2 = 1

# The Chern number of a combined system is the sum of the Chern numbers
# of its constituent parts. This is a fundamental property of topological invariants.
C_total = C1 + C2

# Print the explanation and the calculation.
print("The Chern number of a system formed by joining two subsystems is the sum of their individual Chern numbers.")
print("Given two Chern insulators with Chern numbers C1 and C2:")
print(f"C1 = {C1}")
print(f"C2 = {C2}")
print("The Chern number of the junction (the total system) is C_total = C1 + C2.")
print(f"Therefore, the calculation is:")
print(f"C_total = {C1} + {C2} = {C_total}")