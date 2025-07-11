import math
from fractions import Fraction

# Part 1: Determine n_0 and m_0
# Based on the analysis of heat transfer rates from the plots vs. material thermal conductivities (k_Cu > k_Pt > k_C-steel),
# the lowest heat transfer rate (plot 9) corresponds to Carbon Steel.
n0 = 9

# Based on predicting the heat transfer rate for a carbon steel fin from platinum data,
# the value in plot 9 (7.19 W) is much closer to the prediction for a circular fin (~7.32 W)
# than a square fin (~6.62 W).
# Therefore, the fin is circular.
m0 = 1

# Part 2: Calculate R(c) and R(s)
# For the given conditions, the ratio R = Q_conv / Q_adi simplifies to 1/tanh(mL).

# Calculate R(c) for the circular fin
# Given conditions lead to mL = ln(13)
mL_c = math.log(13)
# R(c) = 1 / tanh(ln(13))
# tanh(ln(13)) = (13^2 - 1) / (13^2 + 1) = 168 / 170 = 84 / 85
Rc = Fraction(1, Fraction(168, 170))

# Calculate R(s) for the square fin
# Given conditions lead to mL = ln(2)
mL_s = math.log(2)
# R(s) = 1 / tanh(ln(2))
# tanh(ln(2)) = (2^2 - 1) / (2^2 + 1) = 3 / 5
Rs = Fraction(1, Fraction(3, 5))

# Part 3: Final Calculation
# The expression to calculate is n0 * (R(c) / R(s))^m0
ratio_R = Rc / Rs
final_result = n0 * (ratio_R)**m0

# Print the final equation with all the numbers
print(f"The calculation is n0 * (R(c) / R(s))^m0")
print(f"Substituting the values:")
print(f"{n0} * (({Rc.numerator}/{Rc.denominator}) / ({Rs.numerator}/{Rs.denominator}))^{m0} = {final_result.numerator}/{final_result.denominator}")