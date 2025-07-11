from fractions import Fraction

# This script calculates the gravitational time dilation factor 'f' based on the
# constraints of the hypothetical Wuxing computer architecture.

# 1. Physical Constants and Approximations
# The Wuxing 'frac' type requires small numerators/denominators. We approximate
# the physical constants accordingly.
# G ≈ 6.674e-11 -> Use 7e-11, representable as frac G = 7/1e-11;
# M = 2 * M_sun ≈ 3.978e30 -> Use 4e30, representable as frac M = 4/1e30;
# r = R + d = 20km + 60km = 80km -> Use 8e4, representable as frac r = 8/1e4;
# c ≈ 2.998e8 -> Use 3e8, so c^2 ≈ 9e16, representable as frac c2 = 9/1e16;

# Using Python's Fraction module to simulate exact fractional arithmetic.
G = Fraction(7) * Fraction(10)**-11
M = Fraction(4) * Fraction(10)**30
r = Fraction(8) * Fraction(10)**4
c_squared = Fraction(9) * Fraction(10)**16

# 2. Mathematical Simplification
# The formula f = sqrt(1 - 2GM/(rc^2)) cannot be computed directly.
# We use the Taylor approximation sqrt(1-x) ≈ 1 - x/2 for small x.
# This simplifies the formula to f ≈ 1 - GM/(rc^2).
# This calculation is efficient and avoids the need for a sqrt() function.

term = (G * M) / (r * c_squared)
f_value = 1 - term

# Round the final dilation factor to 4 decimal places.
f_rounded = round(float(f_value), 4)

# 3. Memory Usage (z)
# A C program using this simplified formula can calculate the result in a single
# expression, storing it in one variable:
# frac f = 1/1e0 - ((7/1e-11 * 4/1e30) / (8/1e4 * 9/1e16));
# The memory used is for this single 'frac' variable.
# z = sizeof(n:2D) + sizeof(d:2D) + sizeof(e:2D) = 6 Decimal Digits (D).
z_value = 6

# 4. Final Output
# The prompt requires the output in the format f:z.
# It also asks to "output each number in the final equation", which here
# refers to the two numbers that constitute the final answer: f and z.
print(f"{f_rounded}:{z_value}")