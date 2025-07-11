import math
import numpy as np

# --- Plan ---
# 1. Define given parameters and constants.
# 2. Determine the subgrade modulus (Es) from the plate bearing test on the subgrade.
# 3. Use the test on the trial section to find the modular ratio (Ep/Es). This requires calculating a deflection factor (F2) and using interpolation on data representing Burmister's charts.
# 4. Calculate the design wheel load force (P_design) and its contact radius (a_design).
# 5. Calculate the required F2 factor (F2_required) for the design load to meet the maximum deflection criteria.
# 6. Use interpolation again with the known modular ratio and F2_required to find the required thickness-to-radius ratio (h/a).
# 7. Calculate the final required pavement thickness (h_design) in mm.
# 8. Print each calculation step clearly showing the formula with the values used.

# --- Step 1: Define given parameters ---
# Plate Bearing Test
P_plate = 30e3  # Load in N (30 kN)
d_plate = 305e-3 # Plate diameter in m (305 mm)
a_plate = d_plate / 2 # Plate radius in m
w_subgrade = 2460e-6 # Deflection on subgrade in m (2460 um)
w_trial = 1080e-6 # Deflection on trial section in m (1080 um)
h_trial = 300e-3 # Trial pavement thickness in m (300 mm)

# Design Parameters
w_allowable = 1.00e-3 # Max allowable deflection in m (1.00 mm)
mass_design = 1.80 * 1000 # Design wheel load mass in kg (1.80 ton)
p_design = 600e3 # Tyre pressure in Pa (600 kN/m2)
g = 9.81 # Acceleration due to gravity in m/s^2

# Material property
mu = 0.5

# Burmister Chart Data for Interpolation
# This data represents the relationship between F2 and Ep/Es for h/a ≈ 2.0
# The F2 values must be in increasing order for numpy.interp
F2_for_ha_2_rev = np.array([0.25, 0.35, 0.44, 0.55, 0.77, 1.0])
Ep_Es_values_rev = np.array([50, 20, 10, 5, 2, 1])

# This data represents the relationship between F2 and h/a for Ep/Es = 10
# The F2 values must be in increasing order for numpy.interp
F2_for_EpEs_10_rev = np.array([0.26, 0.35, 0.44, 0.50, 0.61, 0.81])
ha_ratios_rev = np.array([5.0, 3.0, 2.0, 1.5, 1.0, 0.5])


# --- Step 2: Determine Subgrade Modulus (Es) ---
# For a flexible load with mu=0.5, the surface deflection is w = 1.5 * p * a / E
# where p = P / (pi * a^2)
p_plate = P_plate / (math.pi * a_plate**2)
Es = (1.5 * p_plate * a_plate) / w_subgrade

print("--- Calculation Steps ---")
print("1. Determine Subgrade Modulus (Es)")
print(f"Plate pressure, p_plate = P_plate / (π * a_plate^2) = {P_plate:.2f} N / (π * {a_plate:.4f}^2 m^2) = {p_plate:.2f} Pa")
print(f"Subgrade Modulus, Es = (1.5 * p_plate * a_plate) / w_subgrade")
print(f"Es = (1.5 * {p_plate:.2f} Pa * {a_plate:.4f} m) / {w_subgrade} m = {Es:.2f} Pa ({Es/1e6:.2f} MPa)\n")


# --- Step 3: Determine Modular Ratio (Ep/Es) ---
# For a two-layer system, w = w_subgrade * F2
F2_trial = w_trial / w_subgrade
# Find Ep/Es by interpolating using chart data.
# The ratio h/a for the trial is h_trial / a_plate
h_a_trial = h_trial / a_plate
# Since h_a_trial is very close to 2.0, we use the chart data for h/a ≈ 2.0.
Ep_Es_ratio = np.interp(F2_trial, F2_for_ha_2_rev, Ep_Es_values_rev)
Ep = Es * Ep_Es_ratio

print("2. Determine Pavement Modulus (Ep)")
print(f"Deflection Factor for trial, F2_trial = w_trial / w_subgrade = {w_trial} m / {w_subgrade} m = {F2_trial:.4f}")
print(f"Thickness/Radius ratio for trial, h/a = {h_trial:.4f} m / {a_plate:.4f} m = {h_a_trial:.4f}")
print(f"Using interpolation on Burmister's chart data for h/a ≈ 2, the modulus ratio Ep/Es is found to be {Ep_Es_ratio:.2f}")
print(f"Pavement Modulus, Ep = Es * (Ep/Es) = {Es:.2f} Pa * {Ep_Es_ratio:.2f} = {Ep:.2f} Pa ({Ep/1e6:.2f} MPa)\n")


# --- Step 4: Calculate Design Load Parameters ---
P_design = mass_design * g
# Assuming the tyre pressure is uniform over a circular area: p = P / (pi * a^2)
a_design = math.sqrt(P_design / (math.pi * p_design))

print("3. Determine Design Load Parameters")
print(f"Design Load Force, P_design = mass * g = {mass_design:.2f} kg * {g} m/s^2 = {P_design:.2f} N")
print(f"Contact Radius, a_design = sqrt(P_design / (π * p_design))")
print(f"a_design = sqrt({P_design:.2f} N / (π * {p_design:.2f} Pa)) = {a_design:.4f} m\n")


# --- Step 5: Calculate Required F2 Factor ---
# w_allowable = (1.5 * p_design * a_design / Es) * F2_required
F2_required = (w_allowable * Es) / (1.5 * p_design * a_design)

print("4. Determine Required Deflection Factor (F2_required)")
print(f"F2_required = (w_allowable * Es) / (1.5 * p_design * a_design)")
print(f"F2_required = ({w_allowable} m * {Es:.2f} Pa) / (1.5 * {p_design:.2f} Pa * {a_design:.4f} m) = {F2_required:.4f}\n")


# --- Step 6 & 7: Calculate Required Pavement Thickness ---
# With F2_required and Ep/Es known, find the required h/a ratio by interpolation.
# We assume the curve for Ep/Es=10 is sufficiently close to our calculated Ep/Es=10.11
h_a_ratio_final = np.interp(F2_required, F2_for_EpEs_10_rev, ha_ratios_rev)
h_design = h_a_ratio_final * a_design
h_design_mm = h_design * 1000

print("5. Determine Required Pavement Thickness (h_design)")
print(f"Using interpolation with Ep/Es ≈ 10 and F2_required = {F2_required:.4f}, the required h/a ratio is found to be {h_a_ratio_final:.4f}")
print(f"Required Thickness, h_design = (h/a) * a_design = {h_a_ratio_final:.4f} * {a_design:.4f} m = {h_design:.4f} m")
print("\n--- Final Answer ---")
print(f"The required pavement thickness is {h_design_mm:.2f} mm.")