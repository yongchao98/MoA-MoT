import math
from fractions import Fraction

# Step 1 & 2: Determine n_0 and m_0
# The rate of heat transfer q is proportional to sqrt(k), where k is thermal conductivity.
# Given k_Copper > k_Platinum > k_Carbon_Steel, the q values should be in the same order.
# From the plots: q_1=11.6W, q_7=12.7W (Copper); q_5=8.4W, q_3=7.6W (Platinum); q_9=7.2W (Carbon Steel).
# The lowest q corresponds to plot 9, so n_0 = 9.
# For a given material, a square fin transfers more heat than a circular one.
# For Platinum, q_5 > q_3, so plot 5 is square and plot 3 is circular.
# By comparing the expected heat transfer of a carbon steel fin relative to the platinum fins,
# we find the value in plot 9 (q=7.19W) is consistent with a square cross-section.
# Therefore, the geometry for plot n_0 is square, which means m_0 = -1.

n_0 = 9
m_0 = -1

# Step 3, 4 & 5: Calculate R(c) and R(s)
# The ratio of heat transfer with convection at the tip (Q_conv) to that with an adiabatic tip (Q_adi) is:
# R = Q_conv / Q_adi = (tanh(m*L) + h/(m*k)) / (tanh(m*L) * (1 + h/(m*k)*tanh(m*L)))
# For both the circular and square cases given, the conditions lead to h/(m*k) = 1.
# This simplifies the ratio to R = 1 / tanh(m*L).

# For the circular fin (c):
# The conditions imply m*L = ln(13).
# tanh(ln(13)) = (e^ln(13) - e^-ln(13)) / (e^ln(13) + e^-ln(13))
# tanh(ln(13)) = (13 - 1/13) / (13 + 1/13) = (168/13) / (170/13) = 168/170 = 84/85
Rc_frac = Fraction(85, 84)

# For the square fin (s):
# The conditions imply m*L = ln(2).
# tanh(ln(2)) = (e^ln(2) - e^-ln(2)) / (e^ln(2) + e^-ln(2))
# tanh(ln(2)) = (2 - 1/2) / (2 + 1/2) = (3/2) / (5/2) = 3/5
Rs_frac = Fraction(5, 3)

# Step 6: Final Calculation
# The expression to evaluate is n_0 * (R(c)/R(s))^m_0
ratio_R = Rc_frac / Rs_frac
final_result = n_0 * (ratio_R ** m_0)

print("The equation to be solved is: n_0 * (R(c) / R(s)) ^ m_0")
print(f"Substituting the values, we have:")
print(f"{n_0} * (({Rc_frac.numerator}/{Rc_frac.denominator}) / ({Rs_frac.numerator}/{Rs_frac.denominator})) ^ ({m_0})")
print(f"= {n_0} * ({ratio_R.numerator}/{ratio_R.denominator}) ^ ({m_0})")
print(f"= {n_0} * {ratio_R.denominator}/{ratio_R.numerator}")
print(f"= {n_0 * ratio_R.denominator} / {ratio_R.numerator}")
print(f"The final answer is: {final_result.numerator}/{final_result.denominator}")