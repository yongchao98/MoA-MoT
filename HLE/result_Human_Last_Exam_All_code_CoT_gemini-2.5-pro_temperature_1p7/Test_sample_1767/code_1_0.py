# The Chern number of the first insulator
c1 = 1

# The Chern number of the second insulator
c2 = 1

# The Chern number of a junction of topological insulators is additive.
# Therefore, the total Chern number is the sum of the individual Chern numbers.
c_total = c1 + c2

print(f"The Chern number of the junction is the sum of the individual Chern numbers.")
print(f"Total Chern Number = {c1} + {c2} = {c_total}")