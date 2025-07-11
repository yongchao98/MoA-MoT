import math
from fractions import Fraction

# Step 1: Determine n_0 and m_0 from the problem description and image analysis.
# Based on thermal conductivities (k_Cu > k_Pt > k_CS), the heat transfer rate q
# will follow the same order.
# q values from plots:
# Plot 1: 11.6 W
# Plot 7: 12.7 W
# Plot 5: 8.4 W
# Plot 3: 7.6 W
# Plot 9: 7.2 W
# Grouping by q:
# High q (Copper): 12.7, 11.6
# Mid q (Platinum): 8.4, 7.6
# Low q (Carbon Steel): 7.2
# The carbon steel fin corresponds to the lowest heat transfer, which is in plot 9.
n_0 = 9

# To find the geometry of fin 9, we observe that for both Copper and Platinum,
# one geometry gives ~10% more heat transfer. Let's assume square is higher.
# Q_Pt_s = 8.4, Q_Pt_c = 7.6.
# The ratio of heat transfer Q_CS/Q_Pt should be related to sqrt(k_CS/k_Pt).
# k_CS ~ 50, k_Pt ~ 70. sqrt(50/70) ~ 0.85.
# If fin 9 is square: Q_CS_s ~ 0.85 * Q_Pt_s = 0.85 * 8.4 = 7.14 W. This is close to 7.2 W.
# If fin 9 is circular: Q_CS_c ~ 0.85 * Q_Pt_c = 0.85 * 7.6 = 6.46 W. This is not close.
# Thus, fin 9 is square.
# m_0 = -1 for a square geometry.
m_0 = -1

# Step 2: Calculate R(c) and R(s).
# The ratio R = Q_conv / Q_adi is given by:
# R = (1/tanh(mL)) * (tanh(mL) + h/(mk)) / (1 + (h/(mk))*tanh(mL))
# For both cases given in the problem, it can be shown that h/(mk) = 1.
# When h/(mk) = 1, the formula simplifies to R = 1/tanh(mL).

# For the circular case (c):
# We are given that mL = ln(13).
# tanh(ln(13)) = (e^ln(13) - e^-ln(13)) / (e^ln(13) + e^-ln(13))
# tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = (168/13) / (170/13) = 168/170 = 84/85
tanh_c = Fraction(84, 85)
R_c = 1 / tanh_c

# For the square case (s):
# We are given that mL = ln(2).
# tanh(ln(2)) = (e^ln(2) - e^-ln(2)) / (e^ln(2) + e^-ln(2))
# tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = (3/2) / (5/2) = 3/5
tanh_s = Fraction(3, 5)
R_s = 1 / tanh_s

# Step 3: Compute the final expression.
# The expression is n_0 * (R(c)/R(s))^(m_0)
R_ratio = R_c / R_s
result = n_0 * (R_ratio ** m_0)

# Print the final equation with all numbers substituted
print("The final equation is n_0 * (R(c)/R(s))^m_0")
print(f"Substituting the values:")
print(f"{n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator}))^({m_0})")
print(f"= {n_0} * ({R_ratio.numerator}/{R_ratio.denominator})^({m_0})")
print(f"= {n_0} * ({R_ratio.denominator}/{R_ratio.numerator})")
print(f"= {result.numerator}/{result.denominator}")

# Final answer in the required format
# print(f"\nFinal Answer: {result.numerator}/{result.denominator}")
print(f"\n<<<252/17>>>")