import math
from fractions import Fraction

# Step 1 & 2: Identification of all plots based on physical reasoning.
# The plot indices for g(r) and S(k) are determined by analyzing the shapes
# of the functions, which are dictated by the underlying interaction potentials.

# g(r) indices for the sequence {SS, SR, R, HS, TW}
g_ss = 1  # Plot 1: Repulsive shoulder suppresses g(r) for 1 < r < 1.5
g_sr = 7  # Plot 7: Sticky potential creates a sharp contact peak at r=1
g_r = 3   # Plot 3: Repulsive ramp suppresses g(r) with a piecewise linear shape
g_hs = 9  # Plot 9: Standard hard-sphere g(r) with decaying oscillations
g_tw = 5  # Plot 5: Attractive well enhances g(r) for 1 < r < 1.5

# S(k) indices for the sequence {SS, SR, R, HS, TW}
s_ss = 6  # Plot 6: Repulsive potential leads to low S(0); sharp potential features cause strong oscillations
s_sr = 4  # Plot 4: Strongest attraction (stickiness) gives the highest S(0)
s_r = 8   # Plot 8: Repulsive potential leads to low S(0); smoother potential gives weaker oscillations than SS
s_hs = 0  # S(k) for HS is not plotted (as determined in Step 3)
s_tw = 2  # Plot 2: Attractive potential gives a high S(0)

# Step 3: Identification of the unique system.
# The S(k) for Hard Spheres (HS) is not present among the plots, as none show the expected
# S(0) value for this system (neither the 1D nor 3D value matches any plot).
# Therefore, HS is the unique system, and its S(k) index is 0.

# Step 4: Calculation of R_max for the unique system (HS).
# We use the g(r) for HS from plot 9.
# The formula is R_g(r) = g(r+1)/g(r) for r in {1/2, 3/2, 5/2, ...}.
# The case r=1/2 results in division by zero (g(1/2)=0). We assume the set of r
# for calculation starts at r=3/2.
# We read the values of g(r) from plot 9 at half-integer coordinates (r/σ).
# The grid suggests y-values can be read with reasonable precision.
g_hs_values = {
    1.5: 1.1,
    2.5: 0.9,
    3.5: 1.05,
    4.5: 0.95,
    5.5: 1.02,
}

# Calculate the ratios to find the maximum.
ratios = {
    "R_g(1.5)": g_hs_values[2.5] / g_hs_values[1.5], # ~0.818
    "R_g(2.5)": g_hs_values[3.5] / g_hs_values[2.5], # ~1.167
    "R_g(3.5)": g_hs_values[4.5] / g_hs_values[3.5], # ~0.905
    "R_g(4.5)": g_hs_values[5.5] / g_hs_values[4.5], # ~1.074
}

# The maximum value is R_g(2.5). We express it as a fraction based on the read values.
# R_max = 1.05 / 0.9 = 105 / 90 = 7/6.
R_max_fraction = Fraction(105, 90)
R_max_numerator = R_max_fraction.numerator
R_max_denominator = R_max_fraction.denominator

# Step 5: Format and print the final answer string.
# The sequence is {g(SS), g(SR), g(R), g(HS), g(TW), S(SS), S(SR), S(R), S(HS), S(TW), R_max}
print(f"{{{g_ss},{g_sr},{g_r},{g_hs},{g_tw},"
      f"{s_ss},{s_sr},{s_r},{s_hs},{s_tw},"
      f"{R_max_numerator}/{R_max_denominator}}}")
