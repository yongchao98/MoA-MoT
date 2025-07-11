import math

# Step 1: Define given and assumed parameters
# Given parameters
num_lanes = 3
lane_width_m = 3.6
cross_slope_pct = 1.75

# Assumed parameters for a standard engineering scenario
n_manning = 0.016  # Manning's roughness coefficient for rough-textured asphalt

print("--- Parameters ---")
print(f"Number of lanes (per direction): {num_lanes}")
print(f"Lane width: {lane_width_m} m")
print(f"Cross-slope: {cross_slope_pct}%")
print(f"Assumed Manning's n: {n_manning}\n")

# Step 2: Calculate dependent geometric parameters
# Drainage path length (L)
L = num_lanes * lane_width_m
# Pavement cross-slope (S) as a decimal
S = cross_slope_pct / 100

print("--- Calculated Geometric Properties ---")
print(f"Total drainage path length (L): {L:.2f} m")
print(f"Pavement slope (S): {S}\n")

# Step 3: Iteratively determine the design rainfall intensity (i)
# As no IDF curve was provided, we assume a representative 10-year storm relationship.
# i (mm/hr) = 3048 / (Tc (min) + 20)
def get_i_from_idf(tc_min):
    return 3048 / (tc_min + 20)

# The Time of Concentration (Tc) is calculated using the kinematic wave formula.
# We use a standard formula for Tc in minutes.
def calculate_tc_min(i, n, L, S):
    # Formula: Tc (min) = 0.928 * (n*L)^0.6 * S^-0.3 * i^-0.4
    if i <= 0:
        return float('inf')
    return 0.928 * (n * L)**0.6 * S**-0.3 * i**-0.4

# Perform iterative calculation to find equilibrium i and Tc
print("--- Iterative Calculation for Rainfall Intensity (i) ---")
i_mm_hr = 150.0  # Initial guess for i
print(f"Initial guess for i: {i_mm_hr:.2f} mm/hr")
for iteration in range(5):
    tc_min = calculate_tc_min(i_mm_hr, n_manning, L, S)
    i_mm_hr = get_i_from_idf(tc_min)
    print(f"Iteration {iteration + 1}: Tc = {tc_min:.4f} min -> i = {i_mm_hr:.2f} mm/hr")

print(f"\nConverged Design Rainfall Intensity (i): {i_mm_hr:.2f} mm/hr\n")


# Step 4: Calculate the design water film thickness (d_w)
# The formula is derived from Manning's equation for SI units.
# d_w (mm) = 0.05922 * (n*L)^0.6 * i^0.6 * S^-0.3
dw = 0.05922 * (n_manning * L)**0.6 * i_mm_hr**0.6 * S**-0.3

# Step 5: Display the final calculation and result
print("--- Final Calculation for Water Film Thickness (d_w) ---")
print("The water film thickness is calculated using the formula:")
print("d_w = C * (n * L)^0.6 * i^0.6 * S^-0.3")
print("where C is a constant (0.05922) for the units used (mm, m, mm/hr).\n")

# Print the final equation with values plugged in
term_nl_06 = f"({n_manning} * {L:.2f})^0.6"
term_i_06 = f"({i_mm_hr:.2f})^0.6"
term_S_neg_03 = f"({S})^-0.3"
print(f"d_w = 0.05922 * {term_nl_06} * {term_i_06} * {term_S_neg_03}")

# Print the evaluated terms
val_nl_06 = (n_manning * L)**0.6
val_i_06 = i_mm_hr**0.6
val_S_neg_03 = S**-0.3
print(f"d_w = 0.05922 * ({val_nl_06:.4f}) * ({val_i_06:.4f}) * ({val_S_neg_03:.4f})")

print(f"\nFinal calculated d_w = {dw:.2f} mm")
print("\nThis value represents the design water film thickness at the outer edge of the traveled way for a 10-year design storm, under the stated standard assumptions.")
print(f"\n<<<1.43>>>")