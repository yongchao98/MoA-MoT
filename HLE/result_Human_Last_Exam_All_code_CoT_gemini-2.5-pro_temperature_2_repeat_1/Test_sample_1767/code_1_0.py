# The Chern number is an additive topological quantum number.
# When two Chern insulators are joined together, the Chern number
# of the composite system is the sum of their individual Chern numbers.

# Chern number of the first insulator
c1 = 1

# Chern number of the second insulator
c2 = 1

# Calculate the total Chern number of the junction
c_total = c1 + c2

print(f"The Chern number of the first insulator is {c1}.")
print(f"The Chern number of the second insulator is {c2}.")
print("The total Chern number of the junction is the sum of the two.")
print(f"Final calculation: {c1} + {c2} = {c_total}")