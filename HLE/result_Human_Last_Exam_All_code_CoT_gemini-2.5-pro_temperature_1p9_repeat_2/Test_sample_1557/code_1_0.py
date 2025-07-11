# The task is to identify the correct statements from a list of ten options
# regarding the Duke of Wellington's career.
# Based on historical analysis, the correct statements are:
# 1: His commissariat system from India was adapted for the Peninsular Campaign.
# 6: His methods of integrating local forces became standard colonial practice.
# 8: His use of flying columns was adapted for the Peninsula and later colonial wars.

# The identified correct options are 1, 6, and 8.

# The required output format is a comma-separated string of the
# correct option numbers, sorted in ascending order.

# Define the list of correct option numbers.
correct_options = [1, 6, 8]

# The list is already sorted.
# We will convert each number to a string and then join them with a comma.
# The variable name "final_equation" is used to follow the prompt's instructions,
# even though it's a string of numbers for this particular problem.
final_equation = ",".join(map(str, correct_options))

# Print the final result.
print(final_equation)