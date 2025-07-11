# The Chern number of the first insulator
c1 = 1
# The Chern number of the second insulator
c2 = 1

# The Chern number of a composite system is the sum of the
# Chern numbers of its constituent parts, assuming negligible interaction
# that would alter the bulk topology.
c_total = c1 + c2

print(f"The Chern number of the junction (C_total) is the sum of the individual Chern numbers (C1 and C2).")
print(f"C_total = C1 + C2")
print(f"C_total = {c1} + {c2} = {c_total}")