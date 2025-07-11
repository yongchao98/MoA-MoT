# This script calculates the number of sets (S1, S2, S3, S4) satisfying the given conditions.
# The method is based on a combinatorial argument outlined in the text above.

# Number of choices for placing element 1, based on the condition 1 in S1.
choices_for_1 = 1

# Number of choices for placing element 2, based on the condition 2 in S2.
choices_for_2 = 2

# Number of choices for placing element 3, based on the condition 3 in S3.
choices_for_3 = 3

# Number of choices for placing element 4, for which there are no conditions.
choices_for_4 = 5

# Number of choices for placing element 5, for which there are no conditions.
choices_for_5 = 5

# The total number of ways is the product of the choices for each element.
total_number_of_sets = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

# The final code outputs each number in the final equation.
print(f"The calculation is based on multiplying the number of placement choices for each element:")
print(f"{choices_for_1} (for element 1) * {choices_for_2} (for element 2) * {choices_for_3} (for element 3) * {choices_for_4} (for element 4) * {choices_for_5} (for element 5)")
print(f"Resulting in the final equation:")
print(f"{choices_for_1} * {choices_for_2} * {choices_for_3} * {choices_for_4} * {choices_for_5} = {total_number_of_sets}")