import math

# Step 1: Define the input parameters
# --- Given in the problem ---
lane_width_m = 3.6  # meters
num_lanes = 3
cross_slope_percent = 1.75  # %

# --- Assumed standard design values ---
# Manning's roughness coefficient for rough-textured asphalt
n = 0.016
# Design rainfall intensity for hydroplaning analysis (mm/hr)
i = 200
# Empirical constant for SI units
K = 0.0116

# Step 2: Calculate intermediate values
# L: Total drainage path length in meters
L = num_lanes * lane_width_m
# S: Cross-slope in decimal form (m/m)
S = cross_slope_percent / 100

# Step 3: Calculate the water film thickness using the AASHTO formula
# d_w = [K * (n * L)^0.6 * i^0.4] / S^0.3
term_nL = n * L
term_nL_pow = math.pow(term_nL, 0.6)
term_i_pow = math.pow(i, 0.4)
term_S_pow = math.pow(S, 0.3)

numerator = K * term_nL_pow * term_i_pow
d_w = numerator / term_S_pow

# Step 4: Print the final equation and result
print("Calculation of Water Film Thickness (d_w):\n")
print(f"Formula: d_w = [K * (n * L)^0.6 * i^0.4] / S^0.3\n")
print("Parameters:")
print(f"  K (Constant) = {K}")
print(f"  n (Manning's n) = {n} (Assumed for rough-textured asphalt)")
print(f"  L (Drainage Length) = {L:.1f} m ({num_lanes} lanes * {lane_width_m} m/lane)")
print(f"  i (Rainfall Intensity) = {i} mm/hr (Assumed for hydroplaning design)")
print(f"  S (Cross-slope) = {S} ({cross_slope_percent}%)\n")

print("Final Equation:")
# Displaying the full equation with the numbers plugged in
print(f"d_w (mm) = [{K} * ({n} * {L:.1f})^0.6 * {i}^0.4] / {S}^0.3")
print(f"d_w (mm) = [{K} * ({term_nL:.4f})^0.6 * {i}^0.4] / {S}^0.3")
print(f"d_w (mm) = [{K} * {term_nL_pow:.4f} * {term_i_pow:.4f}] / {term_S_pow:.4f}")
print(f"d_w (mm) = [{numerator:.4f}] / {term_S_pow:.4f}\n")

print(f"Design Water Film Thickness (d_w) = {d_w:.3f} mm")
<<<0.135>>>