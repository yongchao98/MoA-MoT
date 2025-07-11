# Define the Chern numbers for the two insulators
c1 = 1
c2 = 1

# For a composite system made of two Chern insulators with negligible interaction,
# the total Chern number is the sum of the individual Chern numbers.
c_total = c1 + c2

# Print the final equation and the result
print(f"The Chern number of the junction (C_total) is the sum of the individual Chern numbers (C1 and C2).")
print(f"C_total = {c1} + {c2} = {c_total}")