import math

# Step 1: Define the given parameters
L = 40.0
h_target = 10.0

# Step 2: Calculate sin(theta) at the moment when h = 10 m
# From h = L * cos(theta), we have cos(theta) = h / L
cos_theta = h_target / L
# From sin^2(theta) + cos^2(theta) = 1, we find sin(theta)
# The result is positive as the angle is in the first quadrant.
sin_theta = math.sqrt(1 - cos_theta**2)

# Step 3: Use the given value for the rate of change of theta, d(theta)/dt
# The argument to math.cos must be in radians. pi/12 is already in radians.
dtheta_dt = -(3 * math.pi / 10) / math.cos(math.pi / 12)

# Step 4: Calculate the vertical velocity using the formula: dh/dt = -L * sin(theta) * d(theta)/dt
vertical_velocity = -L * sin_theta * dtheta_dt

# Print the breakdown of the calculation as requested
print("The vertical velocity (v) is calculated using the formula: v = -L * sin(theta) * d(theta)/dt")
print("\nThe values used in the final equation are:")
print(f"L = {L}")
print(f"sin(theta) for h=10m is: {sin_theta}")
print(f"d(theta)/dt is given as: {dtheta_dt}")

print("\nSubstituting these numbers into the velocity equation:")
# We use the full variable names here to show each component of the final equation
print(f"v = -({L}) * ({sin_theta}) * ({dtheta_dt})")

# Print the final result
print(f"\nThe calculated vertical velocity of the drawbridge is: {vertical_velocity} m/s")