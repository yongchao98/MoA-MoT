import math
from fractions import Fraction

# Step 1 & 2: Determine n_0 and m_0 from the problem description and image analysis.
# n_0 is the plot number for the carbon steel fin, which has the lowest heat transfer rate.
# The lowest q is in plot 9.
n_0 = 9
# m_0 is 1 for circular geometry, -1 for square.
# Analysis of the image layout suggests plot 9 has a circular cross-section.
m_0 = 1

# Step 3 & 4: Calculate R(c), the ratio for the circular fin.
# The problem provides conditions that simplify the calculation.
# mL = ln(13) and h/(mk) = 1.
# The ratio R = Q_conv / Q_adi simplifies to 1/tanh(mL).
# We use fractions for exact arithmetic.
tanh_ln13 = Fraction(math.sinh(math.log(13)), math.cosh(math.log(13))).limit_denominator()
R_c = 1 / tanh_ln13

# Step 5: Calculate R(s), the ratio for the square fin.
# The problem provides conditions that simplify the calculation.
# mL = ln(2) and h/(mk) = 1.
# The ratio R = Q_conv / Q_adi simplifies to 1/tanh(mL).
tanh_ln2 = Fraction(math.sinh(math.log(2)), math.cosh(math.log(2))).limit_denominator()
R_s = 1 / tanh_ln2

# Step 6: Compute the final expression.
result = n_0 * (R_c / R_s)**m_0

# Print the components of the final equation
print(f"n_0 = {n_0}")
print(f"m_0 = {m_0}")
print(f"R(c) = 1 / tanh(ln(13)) = 1 / ({tanh_ln13.numerator}/{tanh_ln13.denominator}) = {R_c.numerator}/{R_c.denominator}")
print(f"R(s) = 1 / tanh(ln(2)) = 1 / ({tanh_ln2.numerator}/{tanh_ln2.denominator}) = {R_s.numerator}/{R_s.denominator}")
print("\nFinal Calculation:")
print(f"{n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator}))^{m_0} = {result.numerator}/{result.denominator}")