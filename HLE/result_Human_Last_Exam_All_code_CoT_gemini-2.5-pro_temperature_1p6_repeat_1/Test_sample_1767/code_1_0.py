# The Chern number is an additive topological quantum number.
# When two systems are combined, the Chern number of the resulting system is the
# sum of the individual Chern numbers, provided they are weakly coupled.

# Chern number of the first insulator
C1 = 1

# Chern number of the second insulator
C2 = 1

# Calculate the Chern number of the junction
C_total = C1 + C2

# Print the final result as an equation
print(f"The Chern number of the junction is the sum of the individual Chern numbers.")
print(f"C_total = C1 + C2 = {C1} + {C2} = {C_total}")
