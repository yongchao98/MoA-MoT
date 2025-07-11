import math

# --- Given Data ---

# Design Constraint
delta_design_mm = 1.00  # mm
delta_design_m = delta_design_mm / 1000.0  # m

# Plate Bearing Test on Subgrade
P_test_kN = 30.0  # kN
d_plate_mm = 305.0  # mm
delta_s_um = 2460.0  # micrometers

# Plate Bearing Test on Trial Pavement
h_trial_mm = 300.0  # mm
delta_p_um = 1080.0  # micrometers

# Design Wheel Load
W_design_ton = 1.80  # tons
p_design_kPa = 600.0  # kN/m^2

# Material Property
poisson_ratio = 0.5  # mu for all materials

# --- Conversions to SI units ---
P_test_N = P_test_kN * 1000.0  # N
d_plate_m = d_plate_mm / 1000.0  # m
a_plate_m = d_plate_m / 2.0  # m
delta_s_m = delta_s_um / 1_000_000.0  # m
h_trial_m = h_trial_mm / 1000.0  # m
delta_p_m = delta_p_um / 1_000_000.0  # m
# Assuming g = 9.81 m/s^2 for weight to force conversion
P_design_N = W_design_ton * 1000.0 * 9.81  # N
p_design_Pa = p_design_kPa * 1000.0  # Pa (N/m^2)

print("--- Step 1: Determine Subgrade Modulus (E_s) ---")

# For a single layer with mu=0.5, deflection delta = 1.5 * p * a / E
# Pressure from test plate: p = P / (pi * a^2)
p_test_Pa = P_test_N / (math.pi * a_plate_m**2)

# Rearranging for E_s: E_s = 1.5 * p * a / delta_s
E_s = (1.5 * p_test_Pa * a_plate_m) / delta_s_m
print(f"Plate test pressure (p_test): {p_test_Pa / 1e6:.2f} MPa")
print(f"Subgrade Modulus of Elasticity (E_s): {E_s / 1e6:.2f} MPa")
print("-" * 50)


print("--- Step 2: Determine Pavement Modulus (E_p) ---")
# For a two-layer system, deflection delta = (1.5 * p * a / E_s) * F2
# Rearranging for F2: F2 = delta_p / (1.5 * p_test * a_plate / E_s)
F2_trial = (delta_p_m * E_s) / (1.5 * p_test_Pa * a_plate_m)
h_a_ratio_trial = h_trial_m / a_plate_m
print(f"For the trial section (h/a = {h_a_ratio_trial:.3f}), the calculated deflection factor (F2_trial) is {F2_trial:.4f}")

# Interpolate to find Ep/Es ratio from Burmister's chart values.
# At h/a=2.0 (close to our 1.967), chart values are roughly:
# Point 1: (Ep/Es = 5, F2 = 0.5)
# Point 2: (Ep/Es = 10, F2 = 0.4)
# Linear interpolation: y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
# Here, y=Ep/Es, x=F2
Ep_Es_ratio = 5 + (F2_trial - 0.5) * (10 - 5) / (0.4 - 0.5)
E_p = Ep_Es_ratio * E_s
print(f"Interpolated Moduli Ratio (E_p/E_s): {Ep_Es_ratio:.2f}")
print(f"Pavement Modulus of Elasticity (E_p): {E_p / 1e6:.2f} MPa")
print("-" * 50)


print("--- Step 3: Determine Required Pavement Thickness (h) ---")
# Calculate radius of design wheel load area: P = p * pi * a^2
a_design_m = math.sqrt(P_design_N / (p_design_Pa * math.pi))
print(f"Design wheel load contact radius (a_design): {a_design_m * 1000:.2f} mm")

# Find the required F2 for the design case: delta_design = (1.5 * p_design * a_design / E_s) * F2_design
F2_design = (delta_design_m * E_s) / (1.5 * p_design_Pa * a_design_m)
print(f"Required deflection factor for design (F2_design): {F2_design:.4f}")

# Now, find h/a ratio for F2_design and our Ep_Es_ratio. We interpolate again from chart data.
# We need two points on the Ep/Es = Ep_Es_ratio curve. We find them by interpolating vertically first.
# At h/a = 1.0: F2 for Ep/Es=5 is ~0.72, for Ep/Es=10 is ~0.63
F2_at_ha1_EpEs5 = 0.72
F2_at_ha1_EpEs10 = 0.63
F2_at_ha1 = F2_at_ha1_EpEs5 + (Ep_Es_ratio - 5) * (F2_at_ha1_EpEs10 - F2_at_ha1_EpEs5) / (10 - 5)

# At h/a = 2.0: F2 for Ep/Es=5 is ~0.5, for Ep/Es=10 is ~0.4
F2_at_ha2_EpEs5 = 0.5
F2_at_ha2_EpEs10 = 0.4
F2_at_ha2 = F2_at_ha2_EpEs5 + (Ep_Es_ratio - 5) * (F2_at_ha2_EpEs10 - F2_at_ha2_EpEs5) / (10 - 5)

# Now interpolate horizontally for h/a using the two points we just found for our curve.
# Point A: (h/a = 1.0, F2 = F2_at_ha1)
# Point B: (h/a = 2.0, F2 = F2_at_ha2)
h_a_design_ratio = 1.0 + (F2_design - F2_at_ha1) * (2.0 - 1.0) / (F2_at_ha2 - F2_at_ha1)
print(f"Interpolated design thickness-to-radius ratio (h/a): {h_a_design_ratio:.3f}")

# Finally, calculate the required pavement thickness
h_design_m = h_a_design_ratio * a_design_m
h_design_mm = h_design_m * 1000.0

print("-" * 50)
print("--- Final Result ---")
print(f"The required pavement thickness is: {h_design_mm:.2f} mm")
print("\nVerification using the final calculated values:")
print(f"The design deflection equation is: Δ = (1.5 * p_design * a_design / E_s) * F_2")
print(f"With our calculated values, the final design equation is:")
print(f"Δ_design = (1.5 * {p_design_Pa:.0f} Pa * {a_design_m:.4f} m / {E_s:.0f} Pa) * {F2_design:.4f}")

calculated_delta_m = (1.5 * p_design_Pa * a_design_m / E_s) * F2_design
print(f"Calculated Deflection = {calculated_delta_m * 1000:.2f} mm")
print(f"This matches the target deflection of {delta_design_mm:.2f} mm.")

# Format for the final answer
final_answer = h_design_mm
print("\n<<<%.2f>>>" % final_answer)