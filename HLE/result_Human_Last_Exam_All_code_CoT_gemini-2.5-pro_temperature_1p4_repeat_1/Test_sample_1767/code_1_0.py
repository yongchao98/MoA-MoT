# Define the Chern number for the first insulator
C1 = 1

# Define the Chern number for the second insulator
C2 = 1

# When two Chern insulators are joined, the Chern number of the combined system
# is the sum of the individual Chern numbers.
C_total = C1 + C2

# Print the result as an equation
print(f"The Chern number of the junction is the sum of the individual Chern numbers.")
print(f"Total Chern Number = {C1} + {C2} = {C_total}")