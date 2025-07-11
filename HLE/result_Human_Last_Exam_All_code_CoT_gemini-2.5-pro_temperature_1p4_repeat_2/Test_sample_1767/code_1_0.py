# Define the Chern number for each of the two Chern insulators.
chern_number_1 = 1
chern_number_2 = 1

# The Chern number of a composite system is the sum of the Chern numbers
# of its individual parts.
total_chern_number = chern_number_1 + chern_number_2

# Print the explanation and the final equation.
print("The Chern number of the junction formed by two Chern insulators is the sum of their individual Chern numbers.")
print(f"C_total = C1 + C2")
print(f"C_total = {chern_number_1} + {chern_number_2} = {total_chern_number}")