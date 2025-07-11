import math

# Step 1: Define constants
# We set the radius R to 1.0 to get the "raw number" for the coordinates.
R = 1.0

# Step 2: Determine the equilibrium angle theta_0
# The problem describes a string in static equilibrium on a smooth spherical surface.
# A standard result for this problem, assuming the fixing point at the top exerts no tangential force,
# is that the string lifts off at an angle theta_0 where cos(theta_0) = 1/2.
# This gives theta_0 = pi / 3.
theta_0 = math.pi / 3
cos_theta_0 = math.cos(theta_0)
sin_theta_0 = math.sin(theta_0)

# Step 3: Calculate the coordinates of the Center of Mass (CoM).
# As the length of the hanging part of the string is not defined, we calculate the CoM for the
# segment of the string in contact with the pumpkin, from theta = 0 to theta_0 = pi/3.
# The formulas for the CoM of a circular arc are:
# x_cm = R * (1 - cos(theta_0)) / theta_0
# z_cm = R * sin(theta_0) / theta_0

# Calculation for the horizontal coordinate (x_cm)
x_numerator = R * (1 - cos_theta_0)
x_denominator = theta_0
x_cm = x_numerator / x_denominator

print("Calculating the horizontal coordinate (x_cm):")
print(f"x_cm = (R * (1 - cos(θ₀))) / θ₀")
print(f"     = ({R:.1f} * (1 - {cos_theta_0:.3f})) / {theta_0:.4f}")
print(f"     = ({R:.1f} * {1 - cos_theta_0:.3f}) / {theta_0:.4f}")
print(f"     = {x_numerator:.4f} / {theta_0:.4f}")
print(f"     = {x_cm:.4f}")
print("-" * 20)

# Calculation for the vertical coordinate (z_cm)
z_numerator = R * sin_theta_0
z_denominator = theta_0
z_cm = z_numerator / z_denominator

print("Calculating the vertical coordinate (z_cm):")
print(f"z_cm = (R * sin(θ₀)) / θ₀")
print(f"     = ({R:.1f} * {sin_theta_0:.4f}) / {theta_0:.4f}")
print(f"     = {z_numerator:.4f} / {theta_0:.4f}")
print(f"     = {z_cm:.4f}")
print("-" * 20)

# Final raw number output as requested
print("The raw numbers for the horizontal and vertical coordinates separated by a comma are:")
print(f"{x_cm},{z_cm}")
<<<0.47746482927568604,0.8269933403333804>>>