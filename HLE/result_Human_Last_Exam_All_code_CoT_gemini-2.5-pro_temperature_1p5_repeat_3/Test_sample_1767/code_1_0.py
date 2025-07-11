# Define the Chern numbers of the two insulators
c1 = 1
c2 = 1

# The Chern number is an additive property. For a junction of two
# topological materials, the total Chern number is the sum of the
# individual Chern numbers.
c_total = c1 + c2

print("The Chern number of a system composed of two subsystems is the sum of their individual Chern numbers.")
print(f"Chern number of insulator 1 (C1): {c1}")
print(f"Chern number of insulator 2 (C2): {c2}")
print("\nThe Chern number of the junction (C_total) is C1 + C2.")
print(f"The final equation is: {c1} + {c2} = {c_total}")