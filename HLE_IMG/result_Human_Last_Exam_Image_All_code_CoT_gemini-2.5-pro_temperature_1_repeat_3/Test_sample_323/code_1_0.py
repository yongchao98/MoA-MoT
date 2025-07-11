from fractions import Fraction

# Step 1: Determine n_0 and m_0 from the problem description and image analysis.
# Based on the physics of heat transfer, the total heat transfer rate (q) from a fin is primarily
# influenced by the material's thermal conductivity (k). The ranking of thermal conductivities for the
# materials is k_Copper > k_Platinum > k_Carbon_Steel.
# By examining the heat transfer values (q) in plots 1, 3, 5, 7, and 9, we can rank them:
# Plot 7 (q=12.7376 W) > Plot 1 (q=11.634 W) > Plot 5 (q=8.42944 W) > Plot 3 (q=7.62687 W) > Plot 9 (q=7.19061 W).
# The fin with the lowest heat transfer rate must be the one made of carbon steel.
# This corresponds to plot 9.
n_0 = 9

# To determine the geometry (circular or square) of the fin in plot 9, we can compare the ratio of heat transfer
# rates between fins of different materials but the same geometry. This ratio is approximately equal to the
# square root of the ratio of their thermal conductivities: q1/q2 ≈ sqrt(k1/k2).
# Approximate thermal conductivities: k_Platinum ≈ 70 W/mK, k_Carbon_Steel ≈ 50 W/mK.
# Theoretical ratio: sqrt(k_Platinum / k_Carbon_Steel) = sqrt(70/50) ≈ 1.18.
# The fin plots for Platinum are 3 (circular, q=7.62687 W) and 5 (square, q=8.42944 W).
# Let's test the geometry of the fin in plot 9 (Carbon Steel, q=7.19061 W):
# If fin 9 is square (like fin 5): q5/q9 = 8.42944 / 7.19061 ≈ 1.172. This is a very close match to 1.18.
# If fin 9 is circular (like fin 3): q3/q9 = 7.62687 / 7.19061 ≈ 1.061. This is a poor match.
# We conclude that the fin in plot 9 has a square cross-section.
# According to the problem definition, m_0 = -1 for a square geometry.
m_0 = -1

# Step 2: Calculate R(c), the ratio for the circular fin.
# The ratio R is defined as Q_conv / Q_adi.
# For a fin with convection at the tip, Q_conv = M * (tanh(mL) + h/mk) / (1 + (h/mk)tanh(mL)).
# For an adiabatic tip, Q_adi = M * tanh(mL).
# The ratio is R = (tanh(mL) + h/mk) / (tanh(mL) * (1 + (h/mk)tanh(mL))).
# For the circular fin, the given conditions lead to mL = ln(13) and h/mk = 1.
# Substituting h/mk = 1 simplifies the ratio to R = 1 / tanh(mL).
# We calculate tanh(ln(13)) = (13**2 - 1) / (13**2 + 1) = 168 / 170 = 84 / 85.
R_c = Fraction(1) / Fraction(84, 85)

# Step 3: Calculate R(s), the ratio for the square fin.
# The formula for the ratio R remains the same.
# For the square fin, the given conditions lead to mL = ln(2) and h/mk = 1.
# This also simplifies the ratio to R = 1 / tanh(mL).
# We calculate tanh(ln(2)) = (2**2 - 1) / (2**2 + 1) = 3 / 5.
R_s = Fraction(1) / Fraction(3, 5)

# Step 4: Combine the parts to find the final value.
# The expression to evaluate is n_0 * (R(c)/R(s))**m_0.
ratio_R = R_c / R_s
final_value = n_0 * (ratio_R ** m_0)

# Print the final equation with all intermediate values for clarity.
print(f"The final calculation is n_0 * (R(c) / R(s)) ^ m_0")
print(f"= {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator})) ^ ({m_0})")
print(f"= {n_0} * ({ratio_R.numerator}/{ratio_R.denominator}) ^ ({m_0})")
print(f"= {n_0} * ({ratio_R.denominator}/{ratio_R.numerator})")
print(f"= {final_value.numerator}/{final_value.denominator}")
print(f"<<<{final_value.numerator}/{final_value.denominator}>>>")