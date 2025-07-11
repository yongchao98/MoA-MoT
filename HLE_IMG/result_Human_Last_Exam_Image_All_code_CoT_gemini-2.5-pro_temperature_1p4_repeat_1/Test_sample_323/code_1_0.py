from fractions import Fraction

# Step 1: Define n_0 and m_0 based on analysis.
# n_0 is the plot number for the carbon steel fin (lowest heat transfer).
n_0 = 9
# m_0 corresponds to the geometry of the carbon steel fin (square).
m_0 = -1

# Step 2 & 3: Calculate R(c) for the circular fin.
# For the circular case, mL = ln(13) and h/(mk) = 1.
# R(c) simplifies to 1 / tanh(ln(13)).
# tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = 168/170 = 84/85.
R_c = Fraction(85, 84)

# Step 4: Calculate R(s) for the square fin.
# For the square case, mL = ln(2) and h/(mk) = 1.
# R(s) simplifies to 1 / tanh(ln(2)).
# tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = 3/5.
R_s = Fraction(5, 3)

# Step 5: Calculate the final expression.
# The expression is n_0 * (R(c) / R(s))^m_0
ratio_R = R_c / R_s
final_result = n_0 * (ratio_R ** m_0)

# Output the numbers used in the final equation and the result.
print(f"The final expression to calculate is: {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator})) ** ({m_0})")
print(f"The calculated value is: {final_result.numerator}/{final_result.denominator}")
<<<252/17>>>