import math
from fractions import Fraction

# Step 1: Determine n_0, the plot number for the carbon steel fin.
# The heat transfer rate Q is proportional to sqrt(k), where k is thermal conductivity.
# The ranking of thermal conductivities is k_Copper > k_Platinum > k_CarbonSteel.
# Therefore, the fin with the lowest heat transfer rate must be made of carbon steel.
# From the image, the heat transfer rates (q) are:
# q_7 = 12.7376 W (Highest -> Copper)
# q_1 = 11.634 W  (High -> Copper)
# q_5 = 8.42944 W (Mid -> Platinum)
# q_3 = 7.62687 W (Mid -> Platinum)
# q_9 = 7.19061 W (Lowest -> Carbon Steel)
n_0 = 9
print(f"1. Determination of n_0 and m_0:")
print(f"The fin with the lowest heat transfer rate (q = 7.19061 W) is in plot 9. This corresponds to the carbon steel fin.")
print(f"Therefore, n_0 = {n_0}")

# Step 2: Determine m_0, which depends on the geometry of the carbon steel fin.
# We observe that for both copper and platinum, there is a high-Q geometry and a low-Q geometry.
# The ratio Q_high/Q_low is consistent:
# Cu: 12.7376 / 11.634 ~= 1.095
# Pt: 8.42944 / 7.62687 ~= 1.105
# This implies one geometry consistently transfers more heat. Heat transfer Q is proportional to sqrt(Perimeter * Area).
# For comparable dimensions (d~w), a square fin (sqrt(4w * w^2)) transfers more heat than a circular fin (sqrt(pi*d * pi*d^2/4)).
# So, the high-Q geometry is square, and the low-Q geometry is circular.
# To find the geometry of fin 9, we compare it to a platinum fin. The ratio of Q for the same geometry should be sqrt(k_Pt/k_CS).
# k_Pt ~ 71.6, k_CS ~ 54. The theoretical ratio sqrt(71.6/54) is ~1.15.
# Let's check the experimental ratios:
# Ratio with high-Q Pt fin (plot 5, square): 8.42944 / 7.19061 ~= 1.172
# Ratio with low-Q Pt fin (plot 3, circular): 7.62687 / 7.19061 ~= 1.06
# The first ratio is much closer to our theoretical estimate.
# This confirms that fin 9 has the same geometry as fin 5, which is square.
# For a square geometry, m_0 is defined as -1.
m_0 = -1
print(f"The carbon steel fin (plot 9) has a square geometry. Therefore, m_0 = {m_0}\n")

# Step 3: Calculate R(c), the heat transfer ratio for the circular fin.
# The general formula is R = Q_conv/Q_adi = (tanh(mL) + H) / ((1 + H*tanh(mL)) * tanh(mL)), where H = h/(mk).
# For the circular fin, the problem gives conditions that lead to:
# mL = ln(13) and H = h/(mk) = 1.
# When H=1, the formula simplifies to R = 1/tanh(mL).
# tanh(ln(13)) = (e^ln(13) - e^-ln(13)) / (e^ln(13) + e^-ln(13)) = (13 - 1/13) / (13 + 1/13) = 168/170 = 84/85.
Rc = Fraction(85, 84)
print(f"2. Calculation of R(c) and R(s):")
print(f"For the circular fin, R(c) = 1 / tanh(ln(13)) = 1 / (84/85) = {Rc.numerator}/{Rc.denominator}")

# Step 4: Calculate R(s), the heat transfer ratio for the square fin.
# For the square fin, the conditions lead to:
# mL = ln(2) and H = h/(mk) = 1.
# R(s) = 1/tanh(ln(2)).
# tanh(ln(2)) = (e^ln(2) - e^-ln(2)) / (e^ln(2) + e^-ln(2)) = (2 - 1/2) / (2 + 1/2) = 3/5.
Rs = Fraction(5, 3)
print(f"For the square fin, R(s) = 1 / tanh(ln(2)) = 1 / (3/5) = {Rs.numerator}/{Rs.denominator}\n")

# Step 5: Compute the final expression.
# The expression is n_0 * (R(c)/R(s))^m_0
ratio_Rc_Rs = Rc / Rs
final_value = n_0 * (ratio_Rc_Rs) ** m_0

print(f"3. Final Calculation:")
print(f"The expression to evaluate is: n_0 * (R(c) / R(s)) ^ m_0")
# Print the full equation with numbers
print(f"{n_0} * (({Rc.numerator}/{Rc.denominator}) / ({Rs.numerator}/{Rs.denominator})) ^ ({m_0}) = {final_value.numerator}/{final_value.denominator}")
print("\nFinal Answer in x/y format:")
print(f"{final_value.numerator}/{final_value.denominator}")