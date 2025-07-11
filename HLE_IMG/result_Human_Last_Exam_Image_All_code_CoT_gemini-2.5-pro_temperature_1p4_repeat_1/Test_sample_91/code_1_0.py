import math

# Step 1: Define given constants from the problem.
L = 40.0  # Length of the drawbridge in meters
y_target = 10.0  # Target vertical height in meters

# Step 2: Establish the relationship between height y, length L, and angle theta.
# y = L * cos(theta)
# The vertical velocity dy/dt is the time derivative:
# dy/dt = -L * sin(theta) * d(theta)/dt

print("The vertical velocity equation is: dy/dt = -L * sin(theta) * d(theta)/dt")
print(f"Given L = {L} m.")
print("")

# Step 3: Find the values of cos(theta) and sin(theta) when y = 10 m.
print(f"At the moment when y = {y_target} m:")
# From y = L * cos(theta), we find cos(theta)
cos_theta = y_target / L
print(f"cos(theta) = y / L = {y_target} / {L} = {cos_theta}")

# Using sin^2(theta) + cos^2(theta) = 1
sin_theta = math.sqrt(1 - cos_theta**2)
print(f"sin(theta) = sqrt(1 - cos(theta)^2) = sqrt(1 - {cos_theta**2:.4f}) = {sin_theta:.4f}")
print("")

# Step 4: Use the given value for the rate of change of theta, d(theta)/dt.
# d(theta)/dt = -(3 * pi / 10) / cos(pi/12)
pi = math.pi
cos_pi_div_12 = math.cos(pi / 12)
d_theta_dt = -(3 * pi / 10) / cos_pi_div_12

print("The given rate of change for the angle is:")
print(f"d(theta)/dt = -(3 * pi / 10) / cos(pi/12) = {d_theta_dt:.4f} rad/unit_time")
print("")

# Step 5: Substitute all the numerical values into the velocity equation.
print("Substituting the numbers into the velocity equation:")
print(f"dy/dt = -L * sin(theta) * d(theta)/dt")
# The final equation with each number explicitly shown
print(f"dy/dt = -({L}) * ({sin_theta:.4f}) * ({d_theta_dt:.4f})")

# Calculate the final result.
vertical_velocity = -L * sin_theta * d_theta_dt

print(f"\nThe vertical velocity of the drawbridge edge is: {vertical_velocity:.4f} m/unit_time")