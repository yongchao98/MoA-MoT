import math

# The value of the integral is 2**(15/16).
# We calculate this value below.

power_numerator = 15
power_denominator = 16
base = 2

result = base**(power_numerator / power_denominator)

print(f"The value of the integral is {base}^({power_numerator}/{power_denominator})")
print(f"Final calculation: {base}**({power_numerator}/{power_denominator}) = {result}")
