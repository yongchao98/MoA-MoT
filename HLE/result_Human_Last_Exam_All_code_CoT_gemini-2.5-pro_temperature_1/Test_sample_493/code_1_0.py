import math

# This script calculates the average number of stars per constellation based on the derived formula.
# The formula is: Average Size = 8/3 + sqrt(3)/pi

# Define the constants from the formula
pi = math.pi
sqrt3 = math.sqrt(3)

# Calculate the numerator and denominator of the un-simplified fraction: (8*pi + 3*sqrt(3)) / (3*pi)
numerator = 8 * pi + 3 * sqrt3
denominator = 3 * pi

# Calculate the final result
average_constellation_size = numerator / denominator

# As requested, print each number in the final equation
print("The formula for the average constellation size is: (8 * pi + 3 * sqrt(3)) / (3 * pi)")
print(f"Value of pi: {pi}")
print(f"Value of sqrt(3): {sqrt3}")
print("-" * 30)
print(f"Numerator (8 * pi + 3 * sqrt(3)): {numerator}")
print(f"Denominator (3 * pi): {denominator}")
print("-" * 30)
print(f"Average number of stars per constellation = {numerator} / {denominator}")
print(f"Final Answer: {average_constellation_size}")
