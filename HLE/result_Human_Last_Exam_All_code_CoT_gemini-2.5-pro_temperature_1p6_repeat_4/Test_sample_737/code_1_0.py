import math

# --- Given Data ---
# Plate Bearing Test (Rigid Plate)
P_plate = 30 * 1000  # Load in N
d_plate = 305  # Plate diameter in mm
delta_s_um = 2460  # Deflection on subgrade in micrometers
delta_p_um = 1080  # Deflection on trial pavement in micrometers
h_trial = 300  # Trial pavement thickness in mm

# Design Load (Flexible Tire)
W_design_ton = 1.80  # Wheel load in tons
g = 9.81 # acceleration due to gravity in m/s^2
p_design_kPa = 600  # Tyre pressure in kN/m^2
delta_design_max = 1.00  # Maximum allowable deflection in mm

# Material properties
mu = 0.5  # Poisson's ratio for all materials

# --- Calculations ---

# Convert units for consistency (N, mm, MPa)
delta_s = delta_s_um / 1000.0  # mm
delta_p = delta_p_um / 1000.0  # mm
a_plate = d_plate / 2.0  # Plate radius in mm
p_design_MPa = p_design_kPa / 1000.0  # Convert kN/m^2 to MPa (N/mm^2)
P_design = W_design_ton * 1000 * g # Design load in N

print("--- Step 1: Determine Subgrade Modulus (Es) ---")
# For a rigid plate, the average pressure is P / (pi * a^2)
# Deflection formula: delta = 1.18 * p_avg * a / E_s
p_avg_plate = P_plate / (math.pi * a_plate**2)
E_s = (1.18 * p_avg_plate * a_plate) / delta_s
print(f"The formula for subgrade modulus is: E_s = (1.18 * p_avg * a) / delta_s")
print(f"Calculated average pressure on plate (p_avg): {p_avg_plate:.4f} MPa")
print(f"Calculated subgrade modulus (E_s): {E_s:.2f} MPa\n")

print("--- Step 2: Determine Modular Ratio (Ep/Es) ---")
# The deflection factor Fw is the ratio of two-layer deflection to single-layer deflection
Fw_trial = delta_p / delta_s
h_a_ratio_trial = h_trial / a_plate
print(f"Deflection Factor from trial (Fw) = delta_p / delta_s = {delta_p:.4f} / {delta_s:.4f} = {Fw_trial:.4f}")
print(f"h/a ratio for trial section = {h_trial} / {a_plate} = {h_a_ratio_trial:.4f}")

# Simulate Burmister Chart Lookup for Ep/Es
# We will use bilinear interpolation based on standard chart values.
# Known points (Ep/Es, h/a): Fw
# (5, 2): 0.52, (10, 2): 0.42 | (5, 3): 0.42, (10, 3): 0.33
# First, interpolate Fw at h/a=1.9672 for Ep/Es=5 and Ep/Es=10
Fw_at_5_ha2, Fw_at_10_ha2 = 0.52, 0.42
Fw_at_5_ha3, Fw_at_10_ha3 = 0.42, 0.33
Fw_at_Ep5 = Fw_at_5_ha2 + (h_a_ratio_trial - 2) * (Fw_at_5_ha3 - Fw_at_5_ha2) / (3 - 2)
Fw_at_Ep10 = Fw_at_10_ha2 + (h_a_ratio_trial - 2) * (Fw_at_10_ha3 - Fw_at_10_ha2) / (3 - 2)
# Now, interpolate for Ep/Es using the calculated Fw values
Ep_Es_ratio = 5 + (Fw_trial - Fw_at_Ep5) * (10 - 5) / (Fw_at_Ep10 - Fw_at_Ep5)
print(f"By interpolating the Burmister chart, the modular ratio Ep/Es is found to be: {Ep_Es_ratio:.2f}\n")

print("--- Step 3: Determine Required Pavement Thickness (h_design) ---")
# Calculate the radius of the design wheel load contact area
contact_area_design = P_design / p_design_MPa
a_design = math.sqrt(contact_area_design / math.pi)
print(f"Design wheel contact radius (a_design) = sqrt(P_design / (pi * p_design)) = sqrt({P_design:.2f} / (pi * {p_design_MPa:.2f})) = {a_design:.2f} mm")

# Calculate the required deflection factor (Fw) for the design case
# Deflection formula for flexible load: delta = (1.5 * p_design * a_design / E_s) * Fw
deflection_on_subgrade = (1.5 * p_design_MPa * a_design) / E_s
Fw_design = delta_design_max / deflection_on_subgrade
print(f"Deflection on subgrade alone under design load would be: {deflection_on_subgrade:.2f} mm")
print(f"Required Fw to limit deflection to {delta_design_max:.2f} mm = {delta_design_max:.2f} / {deflection_on_subgrade:.2f} = {Fw_design:.4f}\n")

print("--- Step 4: Final Thickness Calculation ---")
# Find h/a ratio from chart for the known Ep/Es ratio and required Fw
# We will interpolate between h/a=3 and h/a=4 using our calculated Ep/Es ratio
Fw_at_ha3_interp = Fw_at_10_ha3 + (Ep_Es_ratio - 10)*(Fw_at_5_ha3 - Fw_at_10_ha3)/(5 - 10)
Fw_at_4_ha_Ep5, Fw_at_4_ha_Ep10 = 0.36, 0.28 # Chart values for h/a = 4
Fw_at_ha4_interp = Fw_at_4_ha_Ep10 + (Ep_Es_ratio - 10)*(Fw_at_4_ha_Ep5 - Fw_at_4_ha_Ep10)/(5-10)

h_a_design = 3 + (Fw_design - Fw_at_ha3_interp) * (4 - 3) / (Fw_at_ha4_interp - Fw_at_ha3_interp)
print(f"Using the chart again with Ep/Es={Ep_Es_ratio:.2f} and Fw={Fw_design:.4f}, the required h/a ratio is: {h_a_design:.2f}")

# Calculate the final thickness
h_design = h_a_design * a_design
print("\nThe required pavement thickness (h_design) is calculated as:")
print(f"h_design = (h/a)_design * a_design = {h_a_design:.2f} * {a_design:.2f}")
print(f"Final Required Pavement Thickness = {h_design:.1f} mm")

print("\n<<<" + f"{h_design:.1f}" + ">>>")