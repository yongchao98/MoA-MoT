# The Chern number is an additive topological invariant.
# When two systems are joined, the Chern number of the combined system
# is the sum of the Chern numbers of the individual systems.

# Chern number of the first insulator
c1 = 1

# Chern number of the second insulator
c2 = 1

# Calculate the Chern number of the junction
c_total = c1 + c2

# Print the calculation and the result
print(f"The Chern number of the junction is the sum of the individual Chern numbers.")
print(f"Total Chern Number = {c1} + {c2} = {c_total}")