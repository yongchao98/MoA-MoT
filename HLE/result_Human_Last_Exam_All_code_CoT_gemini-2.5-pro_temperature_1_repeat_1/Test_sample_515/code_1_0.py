import math

# The point z is (3600, 0). Its squared magnitude is 3600^2.
z_squared_magnitude = 3600.0**2

# The probability is given by the formula derived from potential theory.
# P = 0.5 * (1 - log( (3600^2 + 1) / 3600^2) )
# which simplifies to P = 0.5 * (1 - log(1 + 1/3600^2))

# We express the final equation before calculating the result.
# The calculation shows each number used in the final equation.
print("The final probability is calculated using the equation:")
print("P = 0.5 * (1 - log(1 + 1/3600^2))")
print(f"where log is the natural logarithm, and 3600 is the coordinate from the problem.")
print("")

# Calculate the value
probability = 0.5 * (1 - math.log(1 + 1/z_squared_magnitude))

# Print the result
print(f"The calculated probability is: {probability}")
print(f"The approximate answer with two significant digits is: {probability:.2f}")
