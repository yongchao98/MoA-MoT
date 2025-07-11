import math

# Description:
# This script calculates the horizontal displacement of a guide on a frictionless plane.
# A mass 'm' is released from the top of the guide, slides down a quarter-circle arc,
# across a straight frictional section, and up a second quarter-circle arc.

# --- Given Parameters ---
m = 0.20      # Mass of the body in kg
M = 0.80      # Mass of the guide in kg
R_cm = 20     # Radius of the circular arcs in cm
d_cm = 50     # Length of the straight section in cm
mu_D = 0.20   # Coefficient of dynamic friction

# --- Convert units to SI (meters) ---
R = R_cm / 100.0
d = d_cm / 100.0

print("--- Step 1: Calculate the final height (h) reached by mass m ---")
# Using the work-energy theorem, the loss in potential energy equals the work done by friction.
# m*g*(R-h) = μD*m*g*d  =>  h = R - μD*d
h = R - mu_D * d
print(f"The equation for height h is: h = R - μD * d")
print(f"h = {R:.2f} m - {mu_D:.2f} * {d:.2f} m = {h:.4f} m\n")

print("--- Step 2: Calculate the horizontal displacement of mass m relative to the guide (Δx_g) ---")
# This is the final horizontal coordinate of the mass in the guide's reference frame.
# The initial horizontal coordinate is 0.
# The final coordinate x_g_final is found using the geometry of the second circular arc.
# x_g_final = d + R - sqrt(R^2 - (R-h)^2) which simplifies to d + R - sqrt(2*h*R - h^2)
sqrt_term = math.sqrt(2 * h * R - h**2)
delta_x_g = d + R - sqrt_term
print(f"The equation for the relative displacement Δx_g is: d + R - sqrt(2*h*R - h^2)")
print(f"Δx_g = {d:.2f} + {R:.2f} - sqrt(2*{h:.4f}*{R:.2f} - {h:.4f}^2) = {delta_x_g:.4f} m\n")

print("--- Step 3: Calculate the horizontal displacement of the guide (ΔX_M) ---")
# From the conservation of the center of mass's horizontal position:
# (m + M) * ΔX_M + m * Δx_g = 0
# ΔX_M = - (m / (m + M)) * Δx_g
ratio = m / (m + M)
delta_X_M = -ratio * delta_x_g
print(f"The equation for the guide's displacement ΔX_M is: - (m / (m + M)) * Δx_g")
print(f"ΔX_M = - ({m:.2f} / ({m:.2f} + {M:.2f})) * {delta_x_g:.4f}")
print(f"ΔX_M = -{ratio:.2f} * {delta_x_g:.4f} = {delta_X_M:.4f} m\n")

print("--- Final Answer ---")
print(f"The horizontal displacement of the guide is {delta_X_M:.4f} meters.")
print("The negative sign indicates the guide moves to the left, opposite to the primary horizontal motion of the mass.")

<<< -0.1054 >>>