# Define the Chern number for each insulator
chern_number_1 = 1
chern_number_2 = 1

# The Chern number of the junction is the sum of the individual Chern numbers
total_chern_number = chern_number_1 + chern_number_2

# Print the explanation and the final equation
print("The Chern number of a composite system is the sum of the Chern numbers of its parts.")
print(f"Chern number of the first insulator: C1 = {chern_number_1}")
print(f"Chern number of the second insulator: C2 = {chern_number_2}")
print("The total Chern number of the junction is C_total = C1 + C2.")
print(f"C_total = {chern_number_1} + {chern_number_2} = {total_chern_number}")