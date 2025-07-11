import math

# The problem is to find the furthest distance from point A on the surface of a playdough shape
# optimized to create the strongest gravitational field at A.

# The derived formula for the volume V of this optimal shape is:
# V = (4 * pi * k^3) / 15
# where 'k' is the furthest distance we want to find.

# We are given a volume of 1 cubic meter. We rearrange the formula to solve for k:
# k = (15 * V / (4 * pi))^(1/3)

# The following code calculates the value of k.

# Define the numerical constants from the final equation
volume = 1.0
numerator_constant = 15.0
denominator_constant = 4.0
pi_value = math.pi
power_numerator = 1.0
power_denominator = 3.0

# The calculation is k = (numerator_constant * volume / (denominator_constant * pi)) ^ (power_numerator / power_denominator)
# We print each number used in this equation.
print("To calculate the distance 'k', we use the formula: k = (A / B)^(C)")
print(f"The value for 'A' (15 * V) is: {numerator_constant * volume}")
print(f"The value for 'B' (4 * pi) is: {denominator_constant * pi_value}")
print(f"The value for the exponent 'C' (1/3) is: {power_numerator / power_denominator}")
print("-" * 20)


# Perform the final calculation
result = ((numerator_constant * volume) / (denominator_constant * pi_value)) ** (power_numerator / power_denominator)

print("The distance from point A to the furthest point on the surface of the playdough is (in meters):")
print(result)
