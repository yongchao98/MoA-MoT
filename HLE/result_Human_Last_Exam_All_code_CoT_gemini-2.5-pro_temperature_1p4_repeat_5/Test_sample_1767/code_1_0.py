# The Chern number is an additive topological quantum number.
# When two topological systems, such as Chern insulators, are brought together
# without significant interactions that would close and reopen the bulk band gap,
# the Chern number of the resulting combined system is the sum of the
# Chern numbers of the individual components.

# Define the Chern number for the first insulator
chern_number_1 = 1

# Define the Chern number for the second insulator
chern_number_2 = 1

# The total Chern number of the junction is the sum of the individual Chern numbers.
total_chern_number = chern_number_1 + chern_number_2

print("The Chern number of the first insulator (C1) is: {}".format(chern_number_1))
print("The Chern number of the second insulator (C2) is: {}".format(chern_number_2))
print("\nThe Chern number of the junction is the sum of the individual Chern numbers.")
print("Total Chern Number (C_total) = C1 + C2")

# Display the final calculation and result
print("C_total = {} + {} = {}".format(chern_number_1, chern_number_2, total_chern_number))
