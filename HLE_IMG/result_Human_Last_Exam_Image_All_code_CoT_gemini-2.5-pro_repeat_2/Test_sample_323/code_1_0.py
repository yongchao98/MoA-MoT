import math
from fractions import Fraction

# Step 1: Determine n_0 and m_0 from the problem description and image analysis.
# From analysis, the carbon steel fin corresponds to plot 9.
n_0 = 9
# From the visual layout of the plots, fin 9 is inferred to be circular.
m_0 = 1

# Step 2: Calculate R(c) and R(s) based on the provided formulas and conditions.
# For the given conditions, the ratio R = Q_conv / Q_adi simplifies to 1 / tanh(mL).

# For circular geometry (c), mL = ln(13).
# tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = 168 / 170 = 84 / 85.
R_c = Fraction(85, 84)

# For square geometry (s), mL = ln(2).
# tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = 3/5.
R_s = Fraction(5, 3)

# Step 3: Calculate the final expression.
# The expression to evaluate is n_0 * (R(c)/R(s))^m_0

# Calculate the ratio R(c)/R(s)
ratio_R = R_c / R_s

# Calculate the final value
final_value = n_0 * (ratio_R ** m_0)

# Step 4: Print the components of the equation and the final answer.
print("Calculating the expression: n_0 * (R(c)/R(s))^m_0")
print(f"n_0 = {n_0}")
print(f"m_0 = {m_0}")
print(f"R(c) = {R_c.numerator}/{R_c.denominator}")
print(f"R(s) = {R_s.numerator}/{R_s.denominator}")
print(f"The equation is: {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator}))^{m_0}")
print(f"Result = {final_value.numerator}/{final_value.denominator}")