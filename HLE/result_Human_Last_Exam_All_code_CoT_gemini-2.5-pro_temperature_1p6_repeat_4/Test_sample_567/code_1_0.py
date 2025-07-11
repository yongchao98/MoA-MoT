import math

# This script calculates the value of 'a' at which the only obstruction
# to the symplectic embedding E(1,a) into a 4-ball is the volume constraint.
# This value is known to be a = (7 + 3*sqrt(5)) / 2.

# Define the numbers in the final equation
numerator_constant = 7
sqrt_multiplier = 3
sqrt_operand = 5
denominator = 2

# Perform the calculation
a_value = (numerator_constant + sqrt_multiplier * math.sqrt(sqrt_operand)) / denominator

# As requested, output the final equation with each number and the calculated result.
print(f"The value of a is given by the equation: a = ({numerator_constant} + {sqrt_multiplier} * sqrt({sqrt_operand})) / {denominator}")
print(f"Calculated value: a = {a_value}")
