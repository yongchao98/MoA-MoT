import numpy as np
import math

# --- Step 0: Define constants and given data ---
# Convert all units to a consistent system: Newtons (N) and millimeters (mm)

# Given Poisson's ratio
mu = 0.5

# Plate Bearing Test on Subgrade
P_test = 30 * 1000  # N
d_plate = 305  # mm
a_plate = d_plate / 2  # mm
delta_subgrade = 2460 / 1000  # mm

# Plate Bearing Test on Trial Section
h_trial = 300  # mm
delta_trial = 1080 / 1000  # mm

# Design Wheel Load
mass_design = 1.80 * 1000  # kg
g = 9.81 # m/s^2
P_design = mass_design * g # N
p_tyre_kpa = 600  # kN/m^2
p_tyre_mpa = p_tyre_kpa / 1000 # N/mm^2 (MPa)

# Design Constraint
delta_design_max = 1.00  # mm

print("--- Step 1: Determine Subgrade Modulus (E2) ---")
# Formula for deflection under a rigid plate: Delta = (P * (1 - mu^2)) / (2 * a * E)
# Rearranging for E: E = (P * (1 - mu^2)) / (2 * a * Delta)
E2 = (P_test * (1 - mu**2)) / (2 * a_plate * delta_subgrade)
print(f"The equation for E2 is: (P * (1 - \u03BC\u00b2)) / (2 * a * \u0394)")
print(f"E2 = ({P_test:.0f} N * (1 - {mu:.1f}\u00b2)) / (2 * {a_plate:.1f} mm * {delta_subgrade:.3f} mm)")
print(f"Calculated Subgrade Modulus (E2): {E2:.2f} MPa\n")


print("--- Step 2: Determine Pavement Modulus (E1) ---")
# F2 is the ratio of two-layer deflection to subgrade-only deflection
F2_exp = delta_trial / delta_subgrade
h_a_trial_ratio = h_trial / a_plate
print(f"The trial deflection factor F2 is \u0394_trial / \u0394_subgrade = {delta_trial:.3f} / {delta_subgrade:.3f} = {F2_exp:.4f}")
print(f"The trial thickness/radius ratio (h/a) is {h_trial:.1f} / {a_plate:.1f} = {h_a_trial_ratio:.4f}\n")

# Burmister Chart data for F2 (for rigid plate, mu=0.5)
# This data will be used to interpolate the material properties
h_a_lookup = [0.5, 1.0, 2.0, 4.0, 8.0]
# F2 values for E1/E2 = 5
f2_at_E1_E2_5 = [0.82, 0.67, 0.53, 0.43, 0.38]
# F2 values for E1/E2 = 10
f2_at_E1_E2_10 = [0.73, 0.54, 0.40, 0.31, 0.25]
# F2 values for E1/E2 = 20
f2_at_E1_E2_20 = [0.63, 0.45, 0.31, 0.22, 0.17]

# Find F2 values at h/a = h_a_trial_ratio for different E1/E2 curves
F2_at_h_trial_E5 = np.interp(h_a_trial_ratio, h_a_lookup, f2_at_E1_E2_5)
F2_at_h_trial_E10 = np.interp(h_a_trial_ratio, h_a_lookup, f2_at_E1_E2_10)
F2_at_h_trial_E20 = np.interp(h_a_trial_ratio, h_a_lookup, f2_at_E1_E2_20)

# Interpolate to find the E1/E2 ratio that matches F2_exp
E1_E2_lookup = [5, 10, 20]
F2_lookup_at_h_trial = [F2_at_h_trial_E5, F2_at_h_trial_E10, F2_at_h_trial_E20]
# np.interp requires x-coordinates to be increasing, so we sort by F2
sorted_indices = np.argsort(F2_lookup_at_h_trial)
sorted_F2 = np.array(F2_lookup_at_h_trial)[sorted_indices]
sorted_E1_E2 = np.array(E1_E2_lookup)[sorted_indices]
E1_E2_ratio = np.interp(F2_exp, sorted_F2, sorted_E1_E2)

E1 = E1_E2_ratio * E2
print(f"By interpolating the Burmister chart data for F2 = {F2_exp:.4f} and h/a = {h_a_trial_ratio:.4f}:")
print(f"Calculated Modulus Ratio (E1/E2): {E1_E2_ratio:.2f}")
print(f"Calculated Pavement Modulus (E1) = {E1_E2_ratio:.2f} * {E2:.2f} MPa = {E1:.2f} MPa\n")


print("--- Step 3: Analyze Design Wheel Load ---")
# Find radius 'a' of the design load
contact_area_design = P_design / p_tyre_mpa
a_design = math.sqrt(contact_area_design / math.pi)
print(f"Design wheel load P = {mass_design/1000:.2f} ton * {g} m/s\u00b2 = {P_design:.2f} N")
print(f"Tyre pressure p = {p_tyre_kpa} kN/m\u00b2 = {p_tyre_mpa:.3f} N/mm\u00b2")
print(f"Contact area A = P / p = {P_design:.2f} N / {p_tyre_mpa:.3f} N/mm\u00b2 = {contact_area_design:.2f} mm\u00b2")
print(f"Contact radius a_design = sqrt(A / \u03c0) = {a_design:.2f} mm\n")


print("--- Step 4: Determine Required Pavement Thickness (h) ---")
# Calculate deflection on subgrade-only for the design load
delta_subgrade_design = (P_design * (1 - mu**2)) / (2 * a_design * E2)
# Required F2 factor to meet the max deflection criteria
F2_design = delta_design_max / delta_subgrade_design
print(f"To meet the design criteria, deflection must be \u0394_design = {delta_design_max:.2f} mm")
print(f"The deflection on the subgrade alone under the design load would be {delta_subgrade_design:.3f} mm")
print(f"The required deflection factor F2_design = \u0394_design / \u0394_subgrade_design = {delta_design_max:.2f} / {delta_subgrade_design:.3f} = {F2_design:.4f}\n")

# Now, find the h/a ratio that gives F2_design for our known E1/E2 ratio
# We must first generate the F2 vs h/a curve for our specific E1/E2_ratio
f2_curve_for_design = []
for i in range(len(h_a_lookup)):
    f2_val = np.interp(E1_E2_ratio, E1_E2_lookup, [f2_at_E1_E2_5[i], f2_at_E1_E2_10[i], f2_at_E1_E2_20[i]])
    f2_curve_for_design.append(f2_val)

# Interpolate to find the required h/a ratio
# np.interp requires x-coordinates (F2 values) to be increasing, so we reverse the lists as F2 decreases with h/a
h_a_design_ratio = np.interp(F2_design, f2_curve_for_design[::-1], h_a_lookup[::-1])
h_design = h_a_design_ratio * a_design

print(f"Using the Burmister chart data for E1/E2 = {E1_E2_ratio:.2f} and a required F2 = {F2_design:.4f}:")
print(f"The required thickness/radius ratio (h/a)_design is {h_a_design_ratio:.3f}\n")

print("--- Final Calculation ---")
print(f"The required pavement thickness h is (h/a)_design * a_design.")
print(f"h = {h_a_design_ratio:.3f} * {a_design:.2f} mm")
print(f"h = {h_design:.2f} mm")

print("\n<<<%.2f>>>" % h_design)