# Define the four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# Step 1: Multiply the two 10s
product = num1 * num2

# Step 2: Subtract one 4 from the result
difference = product - num3

# Step 3: Divide the result by the other 4
final_answer = difference / num4

# Print the full equation showing how to get 24
# The int() function is used to display the final answer as a whole number
print(f"(({num1} * {num2}) - {num3}) / {num4} = {int(final_answer)}")