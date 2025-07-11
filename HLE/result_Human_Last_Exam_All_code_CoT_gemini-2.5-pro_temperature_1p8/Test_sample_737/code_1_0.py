import numpy as np

# --- Given Data ---
# Poisson's ratio
mu = 0.5
# Acceleration due to gravity (m/s^2)
g = 9.81

# Plate Bearing Test on Subgrade
P_plate_kN = 30.0  # kN
d_plate_mm = 305.0 # mm
w_s_micro_m = 2460.0 # micrometers

# Plate Bearing Test on Trial Section
h_trial_mm = 300.0 # mm
w_p_micro_m = 1080.0 # micrometers

# Design Wheel Load
W_design_ton = 1.80 # tons
p_design_kPa = 600.0 # kN/m^2
w_design_max_mm = 1.00 # mm

# --- Unit Conversions to SI units (N, m, Pa) ---
P_plate = P_plate_kN * 1000
d_plate = d_plate_mm / 1000
a_plate = d_plate / 2
w_s = w_s_micro_m / 1e6
h_trial = h_trial_mm / 1000
w_p = w_p_micro_m / 1e6
P_design = W_design_ton * 1000 * g
p_design = p_design_kPa * 1000
w_design_max = w_design_max_mm / 1000

# Step 1: Determine the subgrade modulus of elasticity (E_s)
print("--- Step 1: Calculate Subgrade Modulus (E_s) ---")
# Pressure under the plate
A_plate = np.pi * a_plate**2
p_plate = P_plate / A_plate

# Boussinesq's equation for a rigid plate on an elastic half-space
# w = (pi/2) * (1 - mu^2) * (p * a) / E  =>  E = (pi/2) * (1 - mu^2) * (p * a) / w
E_s = (np.pi / 2) * (1 - mu**2) * (p_plate * a_plate) / w_s
print(f"Plate pressure (p_plate): {p_plate/1e6:.2f} MPa")
print(f"Subgrade Modulus of Elasticity (E_s): {E_s/1e6:.2f} MPa\n")


# Step 2: Characterize the two-layer system from the trial section
print("--- Step 2: Characterize Trial Section ---")
# The deflection factor F_w is the ratio of the two-layer deflection to the subgrade-only deflection
F_w_trial = w_p / w_s
h_a_ratio_trial = h_trial / a_plate
print(f"Deflection Factor from Trial (F_w_trial): {F_w_trial:.4f}")
print(f"Thickness/Radius Ratio from Trial (h/a_trial): {h_a_ratio_trial:.4f}")
print("This pair of (h/a, F_w) values defines the pavement/subgrade system for a specific E_p/E_s ratio.\n")

# Step 3: Analyze the design wheel load
print("--- Step 3: Analyze Design Load ---")
# Calculate the contact area and radius for the design load
A_design = P_design / p_design
a_design = np.sqrt(A_design / np.pi)
print(f"Design Load Force (P_design): {P_design:.2f} N")
print(f"Design Contact Radius (a_design): {a_design*1000:.2f} mm\n")

# Step 4: Determine the required pavement thickness
print("--- Step 4: Determine Required Pavement Thickness ---")
# Calculate the theoretical deflection on the subgrade under the DESIGN load
w_s_design = (np.pi / 2) * (1 - mu**2) * (p_design * a_design) / E_s

# Calculate the required deflection factor for the design
F_w_design = w_design_max / w_s_design
print(f"Theoretical subgrade deflection under design load: {w_s_design*1000:.4f} mm")
print(f"Required Deflection Factor for Design (F_w_design): {F_w_design:.4f}")

# Since F_w_design is very close to F_w_trial, we assume the h/a ratio is the same.
# This is a standard engineering approximation in this method.
h_a_ratio_design = h_a_ratio_trial
print(f"Assuming required h/a ratio is same as trial h/a ratio: {h_a_ratio_design:.4f}\n")

# Calculate the final required pavement thickness
h_design = h_a_ratio_design * a_design
h_design_mm = h_design * 1000

print("--- Final Calculation ---")
print("The required pavement thickness (h_design) is calculated by multiplying the design load radius by the h/a ratio determined from the trial section.")
print(f"Required Thickness (mm) = (h/a ratio) * (design radius_m) * 1000")
print(f"{h_design_mm:.2f} mm = {h_a_ratio_design:.4f} * {a_design:.4f} * 1000")
print(f"\nThe required pavement thickness is {h_design_mm:.2f} mm.")

<<<190.44>>>