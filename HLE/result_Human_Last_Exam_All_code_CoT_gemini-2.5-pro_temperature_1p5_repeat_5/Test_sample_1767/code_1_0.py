import sys

# Step 1: Define the Chern numbers for the two insulators.
# According to the problem description, both have a Chern number of 1.
chern_number_1 = 1
chern_number_2 = 1

# Step 2: Calculate the number of chiral edge states at the junction.
# The number of protected chiral edge states at the interface between two
# Chern insulators is given by the difference of their individual Chern numbers.
# This value can be thought of as the effective Chern number of the junction.
junction_chern_number = chern_number_1 - chern_number_2

# Step 3: Print the result in the form of an equation.
# The output will show the calculation step-by-step with the final answer.
print("The effective Chern number of the junction is calculated as the difference between the Chern numbers of the two insulators.")
print(f"C_junction = C1 - C2")
print(f"C_junction = {chern_number_1} - {chern_number_2}")
print(f"C_junction = {junction_chern_number}")

# A result of 0 means the interface is trivial (gapped) and supports no protected chiral edge modes.