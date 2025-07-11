import math

# The value of the integral is 2^(15/16)
base = 2
exponent_num = 15
exponent_den = 16
exponent = exponent_num / exponent_den

result = base ** exponent

# Outputting the equation and the final value as requested
print(f"{base} ** ({exponent_num}/{exponent_den}) = {result}")