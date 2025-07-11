import math
from fractions import Fraction

# Step 1: Determine n_0 and m_0
# Based on thermal conductivity analysis (k_Cu > k_Pt > k_CarbonSteel),
# the fin with the lowest heat transfer for a given geometry is carbon steel.
# Among the circular fins (plots 3, 7, 9), plot 9 has the lowest q (7.19 W).
# So, the carbon steel fin is n_0=9. It's a circular fin, so m_0=1.
n0 = 9
m0 = 1

# Step 2: Calculate R(c) for the circular fin.
# Given mL=ln(13) and h/(mk)=1, R(c) = 1/tanh(ln(13)).
# tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = 168/170 = 84/85.
# R(c) = 1 / (84/85) = 85/84.
Rc_frac = Fraction(85, 84)

# Step 3: Calculate R(s) for the square fin.
# Given mL=ln(2) and h/(mk)=1, R(s) = 1/tanh(ln(2)).
# tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = 3/5.
# R(s) = 1 / (3/5) = 5/3.
Rs_frac = Fraction(5, 3)

# Step 4: Calculate the final expression: n_0 * (R(c)/R(s))^m_0
result_frac = n0 * (Rc_frac / Rs_frac)**m0

# Print the equation with the calculated numbers
print(f"The final calculation is based on the expression: n_0 * (R(c) / R(s))^m_0")
print(f"Substituting the derived values:")
print(f"{n0} * (({Rc_frac.numerator}/{Rc_frac.denominator}) / ({Rs_frac.numerator}/{Rs_frac.denominator}))^{m0}")

# Perform the calculation and show the result as a fraction
print(f"= {n0} * ({Rc_frac / Rs_frac})^{m0}")
print(f"= {n0} * {result_frac / n0}")
print(f"= {result_frac.numerator}/{result_frac.denominator}")