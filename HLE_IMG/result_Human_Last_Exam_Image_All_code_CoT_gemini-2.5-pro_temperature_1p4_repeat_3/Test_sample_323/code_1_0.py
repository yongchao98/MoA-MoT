import math
from fractions import Fraction

# Step 1: Determine n_0 and m_0 from the problem description and image analysis.
# The heat transfer rate Q is proportional to the square root of thermal conductivity (k).
# The order of conductivity is k_copper (~400 W/mK) > k_platinum (~70 W/mK) > k_carbon_steel (~50 W/mK).
# Therefore, we expect Q_copper > Q_platinum > Q_carbon_steel.
# From the plots, the heat transfer rates (Q_conv) are:
# Q_7 = 12.7376 W
# Q_1 = 11.634 W
# Q_5 = 8.42944 W
# Q_3 = 7.62687 W
# Q_9 = 7.19061 W
# The lowest Q value belongs to plot 9, so the carbon steel fin is n_0 = 9.

# To determine its geometry, we observe that for both copper and platinum, there is one fin with higher Q and one with lower Q.
# Plot 7 (Q=12.7) and Plot 1 (Q=11.6) are copper.
# Plot 5 (Q=8.4) and Plot 3 (Q=7.6) are platinum.
# The carbon steel fin (Plot 9, Q=7.19) has a heat transfer rate very close to the square platinum fin (Plot 3, Q=7.6),
# suggesting it shares the same geometry. From visual inspection of the grid, plots with lower Q within a material pair
# (e.g., plot 1 vs 7) are associated with the square fin models.
# Therefore, the carbon steel fin has a square geometry, which means m_0 = -1.

n_0 = 9
m_0 = -1

print("Step 1: Determine parameters from the problem description and image analysis.")
print(f"The carbon steel fin corresponds to plot n_0 = {n_0}.")
print(f"The geometry of this fin is square, so m_0 = {m_0}.")
print("-" * 50)

# Step 2: Calculate the heat transfer ratios R(c) and R(s).
# The ratio of heat transfer with a convective tip to an adiabatic tip is:
# R = Q_conv / Q_adi = (tanh(mL) + h/mk) / (tanh(mL) * (1 + (h/mk)*tanh(mL)))
# For both circular and square cases, the given conditions lead to h/(mk) = 1.
# The formula simplifies to R = 1 / tanh(mL).

print("Step 2: Calculate the heat transfer ratios R(c) and R(s).")

# For the circular fin, R(c):
# The conditions lead to mL = ln(13).
# tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = (168/13) / (170/13) = 168/170 = 84/85
tanh_ln13 = Fraction(84, 85)
Rc = 1 / tanh_ln13
print("For the circular fin:")
print(f"  Given the conditions, mL = ln(13) and h/mk = 1.")
print(f"  tanh(mL) = tanh(ln(13)) = {tanh_ln13.numerator}/{tanh_ln13.denominator}")
print(f"  R(c) = 1 / tanh(mL) = {Rc.numerator}/{Rc.denominator}")
print()

# For the square fin, R(s):
# The conditions lead to mL = ln(2).
# tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = (3/2) / (5/2) = 3/5
tanh_ln2 = Fraction(3, 5)
Rs = 1 / tanh_ln2
print("For the square fin:")
print(f"  Given the conditions, mL = ln(2) and h/mk = 1.")
print(f"  tanh(mL) = tanh(ln(2)) = {tanh_ln2.numerator}/{tanh_ln2.denominator}")
print(f"  R(s) = 1 / tanh(mL) = {Rs.numerator}/{Rs.denominator}")
print("-" * 50)

# Step 3: Calculate the final value.
final_value = n_0 * (Rs / Rc)

print("Step 3: Calculate the final value.")
print(f"The expression to evaluate is: n_0 * (R(c)/R(s))^m_0")
print(f"Substituting the determined values:")
print(f"= {n_0} * (({Rc.numerator}/{Rc.denominator}) / ({Rs.numerator}/{Rs.denominator}))^({m_0})")
print(f"Since m_0 = -1, the expression becomes n_0 * (R(s)/R(c)):")
print(f"= {n_0} * (({Rs.numerator}/{Rs.denominator}) / ({Rc.numerator}/{Rc.denominator}))")
ratio_term = Rs / Rc
print(f"= {n_0} * ({ratio_term.numerator}/{ratio_term.denominator})")
print(f"= ({n_0 * ratio_term.numerator}) / {ratio_term.denominator}")
print(f"The final result is: {final_value.numerator}/{final_value.denominator}")