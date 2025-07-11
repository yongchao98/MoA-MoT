# Number of choices for placing element 1, given 1 in S1
choices_1 = 1

# Number of choices for placing element 2, given 2 in S2
choices_2 = 2

# Number of choices for placing element 3, given 3 in S3
choices_3 = 3

# Number of choices for placing element 4, no constraints
choices_4 = 5

# Number of choices for placing element 5, no constraints
choices_5 = 5

# The total number of ways is the product of the choices for each element
total_number_of_sets = choices_1 * choices_2 * choices_3 * choices_4 * choices_5

# Print the breakdown of the calculation and the final result
print("To find the total number of sets, we multiply the number of choices for each element:")
print(f"Choices for element 1: {choices_1}")
print(f"Choices for element 2: {choices_2}")
print(f"Choices for element 3: {choices_3}")
print(f"Choices for element 4: {choices_4}")
print(f"Choices for element 5: {choices_5}")
print("\nThe final equation is:")
print(f"{choices_1} * {choices_2} * {choices_3} * {choices_4} * {choices_5} = {total_number_of_sets}")
