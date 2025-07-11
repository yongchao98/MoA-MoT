from fractions import Fraction

# Step 1: Determine n_0 and m_0 from the problem description and image analysis.
# The plot with the lowest heat transfer rate corresponds to carbon steel.
# From the image, plot 9 has the lowest q (q=7.19061 W). So, n_0 = 9.
n_0 = 9

# The tip temperature of plot 9's graph (~49 C) matches the temperature in plot 4 (48.6596 C), which shows a square fin.
# Therefore, the geometry is square, and m_0 = -1.
m_0 = -1

# Step 2: Calculate R(c) and R(s).
# The ratio of heat transfer rates is R = Q_conv / Q_adi.
# With the given conditions, it can be shown that the parameter h/(m*k) = 1 for both cases.
# This simplifies the ratio to R = 1 / tanh(mL).

# For the circular fin (c):
# Given: mL = ln(13).
# R(c) = 1 / tanh(ln(13)) = 1 / ((13 - 1/13) / (13 + 1/13)) = 1 / (168/170) = 85/84.
R_c = Fraction(85, 84)

# For the square fin (s):
# Given: mL = ln(2).
# R(s) = 1 / tanh(ln(2)) = 1 / ((2 - 1/2) / (2 + 1/2)) = 1 / (3/5) = 5/3.
R_s = Fraction(5, 3)

# Step 3: Compute the final expression.
# The expression is n_0 * (R(c)/R(s))^m_0.
final_value = n_0 * (R_c / R_s)**m_0

# Print the final equation with all the numbers, as requested.
print(f"The equation to solve is n_0 * (R(c)/R(s))^m_0.")
print("Substituting the derived values:")
print(f"{n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator})) ^ ({m_0}) = {final_value.numerator}/{final_value.denominator}")