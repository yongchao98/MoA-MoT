# Step 1: Define the Chern number for each insulator.
c1 = 1
c2 = 1

# Step 2: Explain the principle.
# The Chern number is an additive topological quantum number.
# When two systems are combined, the total Chern number is the sum
# of the individual Chern numbers.
print("Principle: The Chern number of a composite system is the sum of the Chern numbers of its parts.")
print(f"Chern number of the first insulator (C1): {c1}")
print(f"Chern number of the second insulator (C2): {c2}")
print("-" * 30)

# Step 3: Calculate the total Chern number of the junction.
c_total = c1 + c2

# Step 4: Display the final equation and the result.
# The user requested to see each number in the final equation.
print("The calculation for the total Chern number (C_total) of the junction is:")
print(f"C_total = C1 + C2 = {c1} + {c2}")
print(f"Therefore, the Chern number of the junction is: {c_total}")