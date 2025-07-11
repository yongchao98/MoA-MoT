# The Chern number is an additive topological quantum number.
# When two systems are combined, the Chern number of the resulting system
# is the sum of the individual Chern numbers, provided the junction
# doesn't close the bulk energy gap.

# Chern number of the first insulator
C1 = 1

# Chern number of the second insulator
C2 = 1

# Calculate the Chern number of the junction
C_total = C1 + C2

# Print the result, showing the full equation
print(f"The first Chern insulator has a Chern number C1 = {C1}.")
print(f"The second Chern insulator has a Chern number C2 = {C2}.")
print(f"The Chern number of the junction is the sum of the individual Chern numbers.")
print(f"C_total = C1 + C2")
print(f"{C_total} = {C1} + {C2}")
