# The Chern number is a topological invariant that is additive.
# When two topological systems, like Chern insulators, are brought together
# to form a junction without closing the bulk energy gap, the Chern number
# of the combined system is the sum of the individual Chern numbers.

# Chern number of the first insulator
C1 = 1

# Chern number of the second insulator
C2 = 1

# Calculate the Chern number of the junction
C_total = C1 + C2

print("The Chern number of the junction is the sum of the individual Chern numbers.")
print(f"First insulator's Chern number, C1 = {C1}")
print(f"Second insulator's Chern number, C2 = {C2}")
print("The total Chern number of the junction, C_total, is C1 + C2.")
print(f"C_total = {C1} + {C2} = {C_total}")
