# Define the Chern numbers for the two insulators
chern_insulator_1 = 1
chern_insulator_2 = 1

# The Chern number of the junction is the sum of the individual Chern numbers
# because the Chern number is an additive topological invariant.
chern_junction = chern_insulator_1 + chern_insulator_2

# Print the final equation showing the calculation
print("The Chern number of the first insulator is: 1")
print("The Chern number of the second insulator is: 1")
print(f"The total Chern number of the junction is the sum of the two: {chern_insulator_1} + {chern_insulator_2} = {chern_junction}")