# The Chern number is an additive topological quantum number.
# When two topological systems are joined, the Chern number of the
# combined system is the sum of the individual Chern numbers.

# Chern number of the first insulator
c1 = 1

# Chern number of the second insulator
c2 = 1

# The Chern number of the junction is the sum of the individual Chern numbers.
c_total = c1 + c2

# Print the final equation showing the calculation.
print(f"The total Chern number C_total is C1 + C2.")
print(f"{c1} + {c2} = {c_total}")
