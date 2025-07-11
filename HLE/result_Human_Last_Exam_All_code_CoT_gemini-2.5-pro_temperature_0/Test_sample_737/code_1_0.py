import math

# --- Given Data ---

# Plate Bearing Test on Subgrade
P_plate = 30 * 1000  # Load in N
d_plate = 305  # Plate diameter in mm
a_plate = d_plate / 2  # Plate radius in mm
delta_subgrade = 2460 / 1000  # Deflection on subgrade in mm

# Plate Bearing Test on Trial Pavement
h_trial = 300  # Trial pavement thickness in mm
delta_pavement = 1080 / 1000  # Deflection on trial pavement in mm

# Design Wheel Load
W_design_ton = 1.80  # Weight in tons
g = 9.81 # acceleration due to gravity in m/s^2
P_design = W_design_ton * 1000 * g # Design load in N
q_design_kPa = 600  # Tyre pressure in kN/m^2
q_design_MPa = q_design_kPa / 1000 # Tyre pressure in N/mm^2 (MPa)
delta_design = 1.00  # Maximum allowable deflection in mm

# --- Step 1: Calculate Subgrade Modulus (E_s) ---
print("--- Step 1: Calculate Subgrade Modulus (E_s) ---")
# Formula: E_s = (1.18 * P) / (a * delta)
E_s = (1.18 * P_plate) / (a_plate * delta_subgrade)
print(f"Subgrade test load, P = {P_plate} N")
print(f"Plate radius, a = {a_plate} mm")
print(f"Subgrade deflection, delta_s = {delta_subgrade} mm")
print(f"Calculated Subgrade Modulus, E_s = (1.18 * {P_plate}) / ({a_plate} * {delta_subgrade}) = {E_s:.2f} N/mm^2 (MPa)\n")

# --- Step 2: Determine Pavement to Subgrade Modulus Ratio (E_p/E_s) ---
print("--- Step 2: Determine Pavement to Subgrade Modulus Ratio (E_p/E_s) ---")
# Formula for deflection factor: F2 = delta_pavement / delta_subgrade
F2_trial = delta_pavement / delta_subgrade
# Thickness-to-radius ratio for the trial section
h_a_ratio_trial = h_trial / a_plate
print(f"Trial pavement deflection, delta_p = {delta_pavement} mm")
print(f"Deflection Factor, F2 = {delta_pavement} / {delta_subgrade} = {F2_trial:.4f}")
print(f"Thickness-to-Radius Ratio, h/a = {h_trial} / {a_plate} = {h_a_ratio_trial:.4f}")
# From Burmister's chart, for F2=0.439 and h/a=1.967, the modulus ratio E_p/E_s is approximately 10.
E_p_E_s_ratio = 10
print(f"Based on a standard Burmister chart, the modulus ratio E_p/E_s is determined to be {E_p_E_s_ratio}.\n")

# --- Step 3: Calculate Required Pavement Thickness (h_design) ---
print("--- Step 3: Calculate Required Pavement Thickness (h_design) ---")
# a) Calculate radius of design wheel load
# Formula: a = sqrt(P / (pi * q))
a_design = math.sqrt(P_design / (math.pi * q_design_MPa))
print(f"Design load, P_design = {P_design:.2f} N")
print(f"Design tyre pressure, q_design = {q_design_MPa} N/mm^2")
print(f"Radius of design load, a_design = sqrt({P_design:.2f} / (pi * {q_design_MPa})) = {a_design:.2f} mm\n")

# b) Calculate required deflection factor F2_design
# Formula: F2 = (delta * a * E_s) / (1.18 * P)
F2_design = (delta_design * a_design * E_s) / (1.18 * P_design)
print("b) Calculate the required deflection factor F2_design for the design limit:")
print(f"Required F2_design = ({delta_design} * {a_design:.2f} * {E_s:.2f}) / (1.18 * {P_design:.2f}) = {F2_design:.4f}\n")

# c) Determine required h/a ratio by interpolation from Burmister chart data for E_p/E_s = 10
print("c) Determine the required h/a ratio from the chart data (E_p/E_s = 10):")
# Known points from the chart for E_p/E_s = 10
h_a_1, F2_1 = 1.5, 0.50
h_a_2, F2_2 = 2.0, 0.42
print(f"Using linear interpolation between chart points ({h_a_1}, {F2_1}) and ({h_a_2}, {F2_2}).")
# Linear interpolation formula: y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
# Here, we solve for x (h/a): h_a = h_a_1 + (F2_design - F2_1) * (h_a_2 - h_a_1) / (F2_2 - F2_1)
h_a_ratio_design = h_a_1 + (F2_design - F2_1) * (h_a_2 - h_a_1) / (F2_2 - F2_1)
print(f"Required h/a ratio = {h_a_1} + ({F2_design:.4f} - {F2_1}) * ({h_a_2} - {h_a_1}) / ({F2_2} - {F2_1}) = {h_a_ratio_design:.4f}\n")

# d) Calculate final pavement thickness
print("d) Calculate the final required pavement thickness:")
h_design = h_a_ratio_design * a_design
print(f"Required Pavement Thickness = Required h/a ratio * Radius of design load")
print(f"Required Pavement Thickness = {h_a_ratio_design:.4f} * {a_design:.2f} mm = {h_design:.2f} mm")

# --- Final Answer ---
# The final answer is the required pavement thickness in mm.
final_answer = round(h_design, 2)
print(f"\nThe final required pavement thickness is {final_answer} mm.")