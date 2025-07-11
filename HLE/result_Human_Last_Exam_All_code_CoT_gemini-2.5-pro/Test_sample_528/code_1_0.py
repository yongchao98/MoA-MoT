# Number of choices for placing element 1, given 1 must be in S1.
choices_for_1 = 1

# Number of choices for placing element 2, given 2 must be in S2.
# It can be in S1 or S2\S1.
choices_for_2 = 2

# Number of choices for placing element 3, given 3 must be in S3.
# It can be in S1, S2\S1, or S3\S2.
choices_for_3 = 3

# Number of choices for placing element 4, with no restrictions.
# It can be in any of the 5 regions.
choices_for_4 = 5

# Number of choices for placing element 5, with no restrictions.
# It can be in any of the 5 regions.
choices_for_5 = 5

# The total number of ways is the product of the choices for each element.
total_sets = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

# Print the final equation and the result.
print(f"{choices_for_1} * {choices_for_2} * {choices_for_3} * {choices_for_4} * {choices_for_5} = {total_sets}")