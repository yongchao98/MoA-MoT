import math

# The closed form of the infinite product is 2^(1/8) * exp(-pi/24)

# Define the components of the final equation as requested
base = 2
power_numerator = 1
power_denominator = 8
pi = "pi"
pi_divisor = 24

# Construct the string for the equation
final_equation = f"{base}^({power_numerator}/{power_denominator}) / exp({pi}/{pi_divisor})"

print(f"The closed expression for the infinite product is:\n{final_equation}")

# Calculate the numerical value
numerical_value = math.pow(2, 1/8) * math.exp(-math.pi / 24)

print(f"\nNumerical value: {numerical_value}")