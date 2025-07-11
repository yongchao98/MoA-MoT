import math

# Define the given constants and variables
L = 40.0  # Length of the drawbridge in meters
y = 10.0  # Vertical height of the bridge's edge in meters
pi = math.pi

# Step 1: Find cos(theta) and sin(theta) at the moment y = 10 m
# From the relationship y = L * cos(theta), we have:
# 10 = 40 * cos(theta)
cos_theta = y / L
# Using sin^2(theta) + cos^2(theta) = 1
sin_theta = math.sqrt(1 - cos_theta**2)

# Step 2: Use the given value for the rate of change of theta, d(theta)/dt
# d(theta)/dt = -(3*pi/10) / cos(pi/12)
d_theta_dt = -(3 * pi / 10) / math.cos(pi / 12)

# Step 3: Calculate the vertical velocity (dy/dt) using the formula
# dy/dt = -L * sin(theta) * d(theta)/dt
vertical_velocity = -L * sin_theta * d_theta_dt

# Print the final equation with all the numerical values
print("The vertical velocity (dy/dt) is found using the equation: dy/dt = -L * sin(theta) * d(theta)/dt\n")
print("Substituting the known values:")
print(f"L = {L}")
print(f"sin(theta) = sqrt(1 - (10/40)^2) = {sin_theta:.4f}")
print(f"d(theta)/dt = -(3*pi/10) / cos(pi/12) = {d_theta_dt:.4f} rad/s\n")
print("The final equation with these numbers is:")
print(f"dy/dt = -({L}) * ({sin_theta:.4f}) * ({d_theta_dt:.4f})")

# Print the final result
print(f"\nResult: The vertical velocity of the moving edge is {vertical_velocity:.4f} m/s.")