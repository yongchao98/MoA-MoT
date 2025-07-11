# Define the Chern numbers of the two individual insulators.
c1 = 1
c2 = 1

# When two Chern insulators are joined with a negligible barrier,
# the Chern number of the combined system is the sum of the individual Chern numbers.
c_total = c1 + c2

# Print the explanation and the calculation.
print("The Chern number of the junction is the sum of the individual Chern numbers.")
print(f"Chern number of the first insulator (C1) = {c1}")
print(f"Chern number of the second insulator (C2) = {c2}")
print("The total Chern number (C_total) of the junction is C1 + C2.")
print(f"C_total = {c1} + {c2} = {c_total}")