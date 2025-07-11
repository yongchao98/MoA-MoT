import math
from fractions import Fraction

# Step 1 & 2: Identify n_0 and m_0 from the problem description and image analysis.
# n_0 is the plot number for the carbon steel fin, which has the lowest heat transfer.
# m_0 is based on the geometry of the carbon steel fin (1 for circular, -1 for square).
# Analysis shows plot 9 is the carbon steel fin, and matching its tip temperature with
# the renders confirms it is a square fin.
n_0 = 9
m_0 = -1

# Step 3: Calculate R(c) for the circular geometry condition.
# For the circular case, mL = ln(13). R(c) = 1 / tanh(ln(13)).
# tanh(x) = (e^x - e^-x) / (e^x + e^-x)
# tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = (168/13) / (170/13) = 168/170 = 84/85.
R_c_val_numerator = 84
R_c_val_denominator = 85
R_c = Fraction(R_c_val_denominator, R_c_val_numerator) # 1 / (84/85) = 85/84

# Step 4: Calculate R(s) for the square geometry condition.
# For the square case, mL = ln(2). R(s) = 1 / tanh(ln(2)).
# tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = (3/2) / (5/2) = 3/5.
R_s_val_numerator = 3
R_s_val_denominator = 5
R_s = Fraction(R_s_val_denominator, R_s_val_numerator) # 1 / (3/5) = 5/3

# Step 5: Compute the final expression.
# Expression: n_0 * (R(c) / R(s)) ^ m_0
ratio = R_c / R_s
result = n_0 * (ratio ** m_0)

# Output the equation with the calculated values
print(f"n_0 = {n_0}")
print(f"m_0 = {m_0}")
print(f"R(c) = {R_c.numerator}/{R_c.denominator}")
print(f"R(s) = {R_s.numerator}/{R_s.denominator}")
print(f"Calculation: {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator})) ^ ({m_0})")
print(f"Final result: {result.numerator}/{result.denominator}")
