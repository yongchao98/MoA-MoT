# Number of choices for placing element 1, based on the condition 1 in S_1.
choices_for_1 = 1

# Number of choices for placing element 2, based on the condition 2 in S_2.
# It can be in S_1 or S_2 \ S_1.
choices_for_2 = 2

# Number of choices for placing element 3, based on the condition 3 in S_3.
# It can be in S_1, S_2 \ S_1, or S_3 \ S_2.
choices_for_3 = 3

# Number of choices for placing element 4, which has no constraints.
# It can be in any of the 5 regions: S_1, S_2\S_1, S_3\S_2, S_4\S_3, or U\S_4.
choices_for_4 = 5

# Number of choices for placing element 5, which also has no constraints.
choices_for_5 = 5

# The total number of ways is the product of the choices for each independent element.
total_sets = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

# Print the final equation with each number and the result.
print(f"{choices_for_1} * {choices_for_2} * {choices_for_3} * {choices_for_4} * {choices_for_5} = {total_sets}")