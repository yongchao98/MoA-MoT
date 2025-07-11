import math
from fractions import Fraction

# Step 1: Determine n_0 and m_0 from problem analysis
# Based on the thermal conductivities (Cu > Pt > C-Steel) and geometric factors (Square > Circle),
# we rank the heat transfer rates from the plots:
# Q(Cu,sq) > Q(Cu,circ) > Q(Pt,sq) > Q(Pt,circ) > Q(C-Steel,geom)
# The q values are: 12.7 (P7), 11.6 (P1), 8.4 (P5), 7.6 (P3), 7.2 (P9).
# The lowest heat transfer, 7.19061 W, corresponds to plot 9. Thus, n_0 = 9.
# This value (7.19) is lower than the platinum circular fin's heat transfer (7.63, P3).
# If the carbon steel fin were square, its heat transfer would likely be higher than the platinum circular fin.
# Therefore, the carbon steel fin must have a circular cross-section.
n_0 = 9
# m_0 = 1 for circular geometry.
m_0 = 1

# Step 2: Calculate R(c) and R(s)
# The ratio R is Q_conv / Q_adi = 1 / tanh(mL), when h/(mk) = 1.

# For circular geometry (c):
# Given L*h/k = ln(13) and 4L/d = ln(13).
# We can show that m*L = ln(13) and h/(m*k) = 1.
# R(c) = 1 / tanh(ln(13)).
# tanh(ln(13)) = (e^ln(13) - e^-ln(13)) / (e^ln(13) + e^-ln(13))
# = (13 - 1/13) / (13 + 1/13) = ((169-1)/13) / ((169+1)/13) = 168 / 170 = 84/85.
# R(c) = 1 / (84/85) = 85/84.
Rc = Fraction(85, 84)

# For square geometry (s):
# Given L*h/k = ln(2) and 4L/w = ln(2).
# We can show that m*L = ln(2) and h/(m*k) = 1.
# R(s) = 1 / tanh(ln(2)).
# tanh(ln(2)) = (e^ln(2) - e^-ln(2)) / (e^ln(2) + e^-ln(2))
# = (2 - 1/2) / (2 + 1/2) = (3/2) / (5/2) = 3/5.
# R(s) = 1 / (3/5) = 5/3.
Rs = Fraction(5, 3)

# Step 3: Compute the final value
n0_frac = Fraction(n_0)
ratio_frac = Rc / Rs
result_frac = n0_frac * (ratio_frac ** m_0)

# Print the final equation with all components
print("The final calculation is:")
print(f"{n_0} * (({Rc.numerator}/{Rc.denominator}) / ({Rs.numerator}/{Rs.denominator}))^{m_0} = {result_frac.numerator}/{result_frac.denominator}")