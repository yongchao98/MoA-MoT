# This script identifies and prints the correct options from the historical analysis.

# Based on the evaluation of each statement:
# Statement 1 is correct.
# Statement 6 is correct.
# Statement 8 is correct.

# A list of the correct option numbers.
correct_options = [1, 6, 8]

# The options are already sorted, but sorting is good practice to ensure correctness.
correct_options.sort()

# Convert the list of numbers into a comma-separated string for the final output.
# Each number is included in the final string as requested.
final_answer = ",".join(map(str, correct_options))

print(final_answer)