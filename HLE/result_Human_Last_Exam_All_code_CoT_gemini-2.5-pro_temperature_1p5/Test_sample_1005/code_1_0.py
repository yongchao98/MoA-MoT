import math

# The value of the integral simplifies to 2^(15/16).
# Here we calculate its numerical value.

base = 2
exponent_num = 15
exponent_den = 16
exponent = exponent_num / exponent_den

result = base ** exponent

# Print the result in the format of an equation.
print(f"{base} ** ({exponent_num} / {exponent_den}) = {result}")