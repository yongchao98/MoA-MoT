import math
from fractions import Fraction

# Step 1: Define parameters n_0 and m_0 based on analysis of the provided image.
# n_0 is the plot number for the carbon steel fin. Based on thermal conductivity (k),
# carbon steel has the lowest k, and therefore the lowest heat transfer rate (q).
# q_9 is the lowest value, so n_0 = 9.
n_0 = 9

# m_0 is -1 for a square cross-section and 1 for a circular one.
# Plot 9's tip temperature matches the 3D model in Plot 4, which shows a square fin.
# Thus, m_0 = -1.
m_0 = -1

# Step 2: Calculate R(c) for the circular fin.
# The formula for the ratio R is R = 1 / tanh(mL) when h/(mk) = 1.
# For the circular fin, we are given hL/k = ln(13) and 4L/d = ln(13).
# (mL)^2 = (4L/d)*(hL/k) = ln(13)*ln(13), so mL = ln(13).
# h/(mk) = (hL/k)/mL = ln(13)/ln(13) = 1.
# tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = (168/13) / (170/13) = 168/170 = 84/85.
tanh_mL_c = Fraction(84, 85)
R_c = 1 / tanh_mL_c
# R_c = 85/84

# Step 3: Calculate R(s) for the square fin.
# For the square fin, we are given hL/k = ln(2) and 4L/w = ln(2).
# (mL)^2 = (4L/w)*(hL/k) = ln(2)*ln(2), so mL = ln(2).
# h/(mk) = (hL/k)/mL = ln(2)/ln(2) = 1.
# tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = (3/2) / (5/2) = 3/5.
tanh_mL_s = Fraction(3, 5)
R_s = 1 / tanh_mL_s
# R_s = 5/3

# Step 4: Calculate the final expression: n_0 * (R(c) / R(s)) ^ m_0
# Using the values derived above.
# n_0 = 9, m_0 = -1, R_c = 85/84, R_s = 5/3
# The expression is 9 * ((85/84) / (5/3))^(-1)
# = 9 * ( (85/84) * (3/5) )^(-1)
# = 9 * ( (17*5*3) / (28*3*5) )^(-1)
# = 9 * (17/28)^(-1)
# = 9 * (28/17)
# = 252/17

final_result_numerator = n_0 * 28
final_result_denominator = 17

print(f"The determined parameters are:")
print(f"n_0 = {n_0}")
print(f"m_0 = {m_0}")
print(f"The calculated ratios are:")
print(f"R(c) = {R_c.numerator}/{R_c.denominator}")
print(f"R(s) = {R_s.numerator}/{R_s.denominator}")
print("\nThe final calculation is:")
print(f"value = {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator})) ^ ({m_0})")
print(f"value = {n_0} * ({R_c/R_s}) ^ ({m_0})")
final_value = n_0 * (R_c / R_s)**m_0
print(f"value = {final_value.numerator}/{final_value.denominator}")