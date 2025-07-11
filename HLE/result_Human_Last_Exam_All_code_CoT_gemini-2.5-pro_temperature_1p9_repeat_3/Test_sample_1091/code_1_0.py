import math

# The derived formula for the limit is (2 * sqrt(3)) / pi.
# Here are the numbers used in the final equation.
numerator_factor = 2
sqrt_value = 3
denominator_name = "pi"
denominator_value = math.pi

# Calculate the result
result = (numerator_factor * math.sqrt(sqrt_value)) / denominator_value

# Output the equation and the result.
print(f"The final formula is: ({numerator_factor} * sqrt({sqrt_value})) / {denominator_name}")
print(f"The calculated value of the limit is: {result}")