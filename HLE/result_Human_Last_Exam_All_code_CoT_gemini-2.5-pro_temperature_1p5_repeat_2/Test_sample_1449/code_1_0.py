import math

# Given parameters
m = 100 / 1000  # mass of the block in kg
M = 10.0         # mass of the wedge in kg
h = 2.0          # height of the wedge in meters
theta_deg = 30.0 # angle of the wedge in degrees
mu = 0.5         # coefficient of kinetic friction
g = 10.0         # acceleration due to gravity in m/s^2

# Convert angle to radians for trigonometric functions
theta_rad = math.radians(theta_deg)
sin_theta = math.sin(theta_rad)
cos_theta = math.cos(theta_rad)

# To find the relative acceleration, we first need to find the acceleration of the wedge (A_x).
# This is derived by setting up force equations for both the block and the wedge.
# The formula for the wedge's acceleration (A_x) is:
# A_x = (m * g * cos(theta) * (sin(theta) + mu * cos(theta))) / (M - m * sin(theta)^2 - m * mu * sin(theta) * cos(theta))
numerator_Ax = m * g * cos_theta * (sin_theta + mu * cos_theta)
denominator_Ax = M - m * sin_theta**2 - m * mu * sin_theta * cos_theta
A_x = numerator_Ax / denominator_Ax

# Now, we calculate the acceleration of the block relative to the wedge (a_rel).
# The formula for a_rel, derived from the non-inertial frame analysis, is:
# a_rel = g * (sin(theta) - mu * cos(theta)) - A_x * (cos(theta) + mu * sin(theta))
a_rel = g * (sin_theta - mu * cos_theta) - A_x * (cos_theta + mu * sin_theta)

# Calculate the distance the block slides down the incline
d = h / sin_theta

# Use kinematics to find the time t. d = (1/2) * a_rel * t^2
# t = sqrt(2 * d / a_rel)
time_to_slide = math.sqrt(2 * d / a_rel)

print("Solving for the time it takes the block to slide down the wedge:")
print(f"1. Horizontal acceleration of the wedge (A_x): {A_x:.4f} m/s^2")
print(f"2. Acceleration of the block relative to the wedge (a_rel): {a_rel:.4f} m/s^2")
print(f"3. Distance the block slides along the incline (d): {d:.4f} m")
print("\nFinal Calculation:")
print(f"The time is found using the kinematic equation: t = sqrt(2 * d / a_rel)")
print(f"t = sqrt(2 * {d:.4f} m / {a_rel:.4f} m/s^2)")
print(f"t = {time_to_slide:.4f} s")

# Final answer in the required format
print(f"\n<<<{time_to_slide}>>>")