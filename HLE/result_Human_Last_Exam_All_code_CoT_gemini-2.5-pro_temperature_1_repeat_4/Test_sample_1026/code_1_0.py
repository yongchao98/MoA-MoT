import math

# Define the given parameters
m = 0.20  # Mass of the body in kg
M = 0.80  # Mass of the guide in kg
R_cm = 20  # Radius of the circular arcs in cm
d_cm = 50  # Length of the straight section in cm
mu_D = 0.20  # Coefficient of friction

# Convert units from cm to meters
R = R_cm / 100.0
d = d_cm / 100.0

# --- Step 1: Calculate the final height h ---
# The final height h is determined by the energy lost to friction.
# The work-energy equation is: m*g*h - m*g*R = -mu_D*m*g*d
# Simplifying gives: h = R - mu_D * d
h = R - mu_D * d

# --- Step 2: Calculate the relative horizontal displacement of the mass ---
# This is the horizontal distance the mass travels relative to the guide.
# Initial x-position (relative to guide's start): x_i = 0
# Final x-position (relative to guide's start): x_f = R + d - sqrt(R^2 - (R-h)^2)
# The displacement is delta_x_m_rel = x_f - x_i
delta_x_m_rel = R + d - math.sqrt(R**2 - (R - h)**2)

# --- Step 3: Calculate the horizontal displacement of the guide ---
# From conservation of the center of mass position, we have:
# delta_x_M = - (m / (m + M)) * delta_x_m_rel
delta_x_M = - (m / (m + M)) * delta_x_m_rel

# --- Output the results with numbers in the equations ---
print("This script calculates the horizontal displacement of the guide.")
print("\n--- Calculation Breakdown ---")

print("1. Final height (h) reached by mass m:")
print(f"h = R - mu_D * d")
print(f"h = {R:.2f} m - {mu_D:.2f} * {d:.2f} m = {h:.2f} m")

print("\n2. Horizontal displacement of the mass relative to the guide (Δx_rel):")
print(f"Δx_rel = R + d - sqrt(R² - (R-h)²)")
print(f"Δx_rel = {R:.2f} + {d:.2f} - sqrt({R:.2f}² - ({R:.2f}-{h:.2f})²) = {delta_x_m_rel:.4f} m")

print("\n3. Final horizontal displacement of the guide (Δx_M):")
print(f"Δx_M = - (m / (m + M)) * Δx_rel")
print(f"Δx_M = - ({m:.2f} kg / ({m:.2f} kg + {M:.2f} kg)) * {delta_x_m_rel:.4f} m")
print(f"Δx_M = {delta_x_M:.4f} m")

# Convert final answer to cm for better readability
final_result_cm = delta_x_M * 100
print(f"\nThe horizontal displacement of the guide is {delta_x_M:.4f} meters, or {final_result_cm:.2f} cm.")
print("The negative sign indicates the guide moves to the left.")

# Final answer in the required format (in meters)
print(f"\n<<<{delta_x_M:.4f}>>>")