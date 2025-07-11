import math
from fractions import Fraction

# Step 1 & 2: Determine n_0 and m_0 from the problem description and plot analysis.
# Based on the reasoning that the Carbon Steel fin has the lowest heat transfer rate
# among fins of similar geometry and that it's identifiable as plot 9.
# Further analysis of heat transfer rates and the R=Q_conv/Q_adi ratio
# reveals that fin 9 has a square cross-section.
# Therefore, n_0 = 9 and m_0 = -1.
n_0 = 9
m_0 = -1

# Step 3: Calculate R(c) and R(s)
# The general formula for the ratio R = Q_conv / Q_adi is:
# R = (1 + H/t) / (1 + H*t), where t = tanh(mL) and H = h/(mk)

# For the circular fin (c):
# Given: hL/k = ln(13) and 4L/d = ln(13).
# We can derive mL_c = ln(13) and H_c = 1.
t_c_val = math.tanh(math.log(13))
# For exact fraction, t_c = (e^ln(13) - e^-ln(13)) / (e^ln(13) + e^-ln(13)) = (13 - 1/13) / (13 + 1/13) = (168/13)/(170/13) = 84/85
t_c = Fraction(84, 85)
H_c = 1
R_c = (1 + H_c / t_c) / (1 + H_c * t_c)

# For the square fin (s):
# Given: hL/k = ln(2) and 4L/w = ln(2).
# We can derive mL_s = ln(2) and H_s = 1.
t_s_val = math.tanh(math.log(2))
# For exact fraction, t_s = (e^ln(2) - e^-ln(2)) / (e^ln(2) + e^-ln(2)) = (2 - 1/2) / (2 + 1/2) = (3/2)/(5/2) = 3/5
t_s = Fraction(3, 5)
H_s = 1
R_s = (1 + H_s / t_s) / (1 + H_s * t_s)

# Step 4: Calculate the final value: n_0 * (R(c)/R(s))^m_0
base = R_c / R_s
final_value = n_0 * (base ** m_0)

# Output the components of the final equation and the result.
print(f"The equation to solve is: n_0 * (R(c)/R(s))^m_0")
print(f"From analysis, we found n_0 = {n_0}")
print(f"The geometry for fin n_0 is square, so m_0 = {m_0}")
print(f"The calculated value for R(c) is {R_c.numerator}/{R_c.denominator}")
print(f"The calculated value for R(s) is {R_s.numerator}/{R_s.denominator}")
print(f"The final calculation is: {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator}))^({m_0})")
print(f"The result is: {final_value.numerator}/{final_value.denominator}")