import math

# --- Given Data ---

# Plate Bearing Test
P_plate = 30e3  # Load on plate (N)
d_plate = 305  # Diameter of plate (mm)
a_plate = (d_plate / 2) / 1000 # Radius of plate (m)
delta_subgrade = 2460e-6 # Deflection on subgrade (m)
delta_trial = 1080e-6  # Deflection on trial section (m)
h_trial = 300 / 1000 # Thickness of trial pavement (m)

# Design Load
W_design_ton = 1.80 # Design wheel load (ton)
P_design = W_design_ton * 1000 * 9.81 # Convert mass (ton) to force (N)
p_design = 600e3  # Tyre pressure (N/m^2)

# Design Limit
delta_design_limit = 1.00 / 1000 # Maximum allowed deflection (m)

# --- Step A: Determine Subgrade Modulus (E2) ---

# For a single layer (subgrade test), the deflection is Δ = 1.5 * P / (π * E₂ * a)
# Rearranging for E₂: E₂ = 1.5 * P / (π * Δ * a)
E2 = (1.5 * P_plate) / (math.pi * delta_subgrade * a_plate)

print("--- Step A: Subgrade Properties ---")
print(f"Subgrade Modulus of Elasticity (E₂) = {E2 / 1e6:.2f} MPa\n")

# --- Step B: Characterize the Pavement System from Trial Section ---

# The deflection factor F₂ is the ratio of two-layer deflection to single-layer deflection.
# F₂ = Δ_trial / Δ_subgrade
F2_trial = delta_trial / delta_subgrade
h_over_a_trial = h_trial / a_plate

print("--- Step B: Trial Section Analysis ---")
print(f"The trial pavement gives a deflection factor (F₂_trial) of {F2_trial:.4f}")
print(f"at a thickness/radius ratio (h/a_trial) of {h_over_a_trial:.4f}.")
print("This pair of values corresponds to a specific material stiffness ratio (E₁/E₂) on the Burmister chart.\n")

# --- Step C: Analyze the Design Load ---

# The contact area is assumed to be a circle: Area = P_design / p_design
contact_area_design = P_design / p_design
# Area = π * a² => a = sqrt(Area / π)
a_design = math.sqrt(contact_area_design / math.pi)

print("--- Step C: Design Load Parameters ---")
print(f"Design wheel load (P_design) = {P_design / 1000:.2f} kN")
print(f"Radius of design load contact area (a_design) = {a_design * 1000:.2f} mm\n")

# --- Step D: Determine Required Pavement Thickness (h_design) ---

# We use the main deflection formula to find the required F₂ for the design case.
# Δ = 1.5 * (p * a / E₂) * F₂
# Rearranging for F₂: F₂_design = (Δ_design * E₂) / (1.5 * p_design * a_design)
F2_design = (delta_design_limit * E2) / (1.5 * p_design * a_design)

print("--- Step D: Required Pavement Thickness Calculation ---")
print(f"The required deflection factor (F₂_design) to limit deflection to {delta_design_limit*1000:.2f} mm is {F2_design:.4f}.")

# Compare F₂_design with F₂_trial. Since they are nearly identical, the required h/a ratio
# for the design case will be nearly identical to the h/a ratio from the trial.
print(f"Comparing F₂_design ({F2_design:.4f}) with F₂_trial ({F2_trial:.4f}), we see they are almost identical.")
print("Therefore, we can approximate that the required h/a ratio is the same as the trial's ratio.\n")
h_over_a_design = h_over_a_trial # This is our key approximation

# Calculate the final thickness h_design = (h/a)_design * a_design
h_design = h_over_a_design * a_design

print("--- Final Answer ---")
print("The required pavement thickness (h_design) is calculated as:")
print(f"h_design (mm) = (h/a ratio for design) * a_design (mm)")
print(f"h_design (mm) = {h_over_a_design:.4f} * {a_design * 1000:.2f}")
print(f"h_design (mm) = {h_design * 1000:.2f}")

final_thickness_mm = h_design * 1000
# The final result in the required format.
print(f"\n<<<{final_thickness_mm:.2f}>>>")