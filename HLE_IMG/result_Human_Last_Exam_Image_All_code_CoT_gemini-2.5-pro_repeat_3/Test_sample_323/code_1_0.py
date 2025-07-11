import math
from fractions import Fraction

# Step 1 & 2: Determine n_0 and m_0
# From the problem description, we know k_Copper > k_Platinum > k_Carbon Steel.
# The heat transfer rate 'q' is proportional to some function of k. Higher k generally leads to higher q.
# The heat transfer values (q_conv) from the plots are:
# Plot 1: 11.634 W
# Plot 3: 7.62687 W
# Plot 5: 8.42944 W
# Plot 7: 12.7376 W
# Plot 9: 7.19061 W
# Ordering these: q7 > q1 > q5 > q3 > q9.
# The lowest heat transfer rate is in plot 9. Therefore, plot 9 represents the carbon steel fin.
n_0 = 9
print(f"Analysis of heat transfer rates indicates the carbon steel fin is in plot n_0 = {n_0}.")

# To determine the geometry of fin 9 (m_0), we observe the shape of the temperature curves.
# By pairing the data plots with visualization plots (based on tip temperatures), we find:
# Plots 1 & 3 correspond to square fins.
# Plots 5 & 7 correspond to circular fins.
# The curve shape of plot 9 is visually most similar to the curves of plots 1 and 3.
# We therefore deduce that the fin in plot 9 has a square cross-section.
# For a square geometry, m_0 is defined as -1.
m_0 = -1
print(f"Visual analysis of the curve shape suggests the fin has a square geometry, so m_0 = {m_0}.")
print("-" * 20)

# Step 3: Calculate R(c) and R(s)
# The ratio R is given by R = (Q_conv / Q_adi) = 1 / tanh(mL) under the specific conditions
# provided in the problem, as for both cases h/(mk) simplifies to 1.

# For the circular fin (c):
# Given mL = ln(13).
# R(c) = 1 / tanh(ln(13))
# tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = (168/13) / (170/13) = 168/170 = 84/85.
# So, R(c) = 85/84.
Rc_num = 85
Rc_den = 84
Rc = Fraction(Rc_num, Rc_den)
print(f"The ratio for the circular fin is R(c) = 1 / tanh(ln(13)) = {Rc.denominator}/{Rc.numerator} = {Rc_num}/{Rc_den}.")

# For the square fin (s):
# Given mL = ln(2).
# R(s) = 1 / tanh(ln(2))
# tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = (3/2) / (5/2) = 3/5.
# So, R(s) = 5/3.
Rs_num = 5
Rs_den = 3
Rs = Fraction(Rs_num, Rs_den)
print(f"The ratio for the square fin is R(s) = 1 / tanh(ln(2)) = {Rs.denominator}/{Rs.numerator} = {Rs_num}/{Rs_den}.")
print("-" * 20)

# Step 4: Calculate the final expression
# The expression is n_0 * (R(c)/R(s))^m_0
ratio_Rc_Rs = Rc / Rs
final_value = n_0 * (ratio_Rc_Rs) ** m_0

print("Final Calculation:")
print(f"The expression to evaluate is: n_0 * (R(c) / R(s)) ^ m_0")
print(f"= {n_0} * (({Rc_num}/{Rc_den}) / ({Rs_num}/{Rs_den})) ^ ({m_0})")
print(f"= {n_0} * ({ratio_Rc_Rs.numerator}/{ratio_Rc_Rs.denominator}) ^ ({m_0})")
if m_0 == -1:
    inv_ratio = Fraction(ratio_Rc_Rs.denominator, ratio_Rc_Rs.numerator)
    print(f"= {n_0} * ({inv_ratio.numerator}/{inv_ratio.denominator})")
    print(f"= ({n_0 * inv_ratio.numerator}) / {inv_ratio.denominator}")

print(f"= {final_value.numerator}/{final_value.denominator}")
