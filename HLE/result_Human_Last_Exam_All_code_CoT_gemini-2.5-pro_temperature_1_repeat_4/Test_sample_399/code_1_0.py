from fractions import Fraction

# After solving the differential equation with the initial condition,
# we get the particular solution:
# x(t) = 1 / (-t^2/2 + t - 1/8)
# We need to find the value of x(1).

# Define the numbers that appear in the final calculation for x(1):
# x(1) = num_1 / (-(num_2**num_3 / num_4) + num_5 - (num_6 / num_7))
num_1 = 1  # Numerator
num_2 = 1  # Value of t
num_3 = 2  # Power of t
num_4 = 2  # Denominator for the t term
num_5 = 1  # The additive t term
num_6 = 1  # Numerator of the constant C
num_7 = 8  # Denominator of the constant C

# Calculate the result using the Fraction module for an exact answer.
denominator = -(Fraction(num_2**num_3, num_4)) + Fraction(num_5) - Fraction(num_6, num_7)
result = Fraction(num_1, denominator)

# Print the final equation with all its components, and the final answer.
print("The equation to calculate x(1) is derived from the particular solution evaluated at t=1.")
print(f"x(1) = {num_1} / (-({num_2}^{num_3}/{num_4}) + {num_5} - ({num_6}/{num_7}))")
print(f"x(1) = 1 / (-1/2 + 1 - 1/8)")
print(f"x(1) = 1 / (1/2 - 1/8)")
print(f"x(1) = 1 / (3/8)")
print(f"The final value of x(1) is: {result}")