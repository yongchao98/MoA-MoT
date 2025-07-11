import math

# Given values from the problem
L = 40.0  # Length of the drawbridge in meters
y = 10.0  # Vertical height of the edge in meters

# Step 1: Find cos(theta) and sin(theta) at the given instant.
# The geometric relationship is y = L * cos(theta).
cos_theta = y / L
# Using the identity sin^2(theta) + cos^2(theta) = 1, we find sin(theta).
# We take the positive root because the angle theta is acute.
sin_theta = math.sqrt(1 - cos_theta**2)

# Step 2: Calculate the given rate of change of theta, d(theta)/dt.
cos_pi_12 = math.cos(math.pi / 12)
dtheta_dt = (-3 * math.pi) / (10 * cos_pi_12)

# Step 3: Calculate the vertical velocity, dy/dt.
# The velocity relationship is dy/dt = -L * sin(theta) * d(theta)/dt.
dy_dt = -L * sin_theta * dtheta_dt

# Print the process and the final equation with values substituted.
print("The relationship between the vertical height (y) and the angle (theta) is: y = L * cos(theta)")
print("The vertical velocity (dy/dt) is found by differentiating with respect to time: dy/dt = -L * sin(theta) * d(theta)/dt")
print("\nAt the moment the bridge is 10m above the ground:")
print(f"L = {L}")
print(f"y = {y}")
print(f"cos(theta) = y / L = {y} / {L} = {cos_theta}")
print(f"sin(theta) = sqrt(1 - cos(theta)^2) = sqrt(1 - {cos_theta:.2f}^2) = {sin_theta:.4f}")
print(f"The given rate of change of the angle is d(theta)/dt = -3*pi / (10 * cos(pi/12)) = {dtheta_dt:.4f}")
print("\nSubstituting these into the velocity equation:")
print(f"dy/dt = -({L}) * ({sin_theta:.4f}) * ({dtheta_dt:.4f})")
print(f"\nThe vertical velocity of the moving edge is {dy_dt:.4f} m/unit_time.")

print(f"\nFinal calculated value: {dy_dt}")
<<<37.78843513364966>>>