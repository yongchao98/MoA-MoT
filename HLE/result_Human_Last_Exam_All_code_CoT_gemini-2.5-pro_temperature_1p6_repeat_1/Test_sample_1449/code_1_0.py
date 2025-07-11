import math

# Given parameters
m = 0.1  # kg (mass of the block)
M = 10.0 # kg (mass of the wedge)
theta_deg = 30.0 # degrees (angle of the wedge)
h = 2.0  # m (initial height)
mu = 0.5 # coefficient of friction
g = 10.0 # m/s^2 (acceleration due to gravity)

# Convert angle to radians
theta_rad = math.radians(theta_deg)
sin_theta = math.sin(theta_rad)
cos_theta = math.cos(theta_rad)

# Step 1: Calculate the acceleration of the wedge (A)
# Using Newton's second law for both the block and the wedge, we can derive an expression for A.
# The horizontal force on the wedge from the block causes it to accelerate.
# Solving the system of equations from the force analysis yields:
numerator_A = m * g * cos_theta * (sin_theta + mu * cos_theta)
denominator_A = M + m * sin_theta**2 + m * mu * sin_theta * cos_theta
A = numerator_A / denominator_A

# Step 2: Calculate the relative acceleration of the block down the wedge (a)
# This is derived from the forces on the block in the non-inertial frame of the wedge.
# a = g*(sin(theta) - mu*cos(theta)) - A*(cos(theta) - mu*sin(theta))
a = g * (sin_theta - mu * cos_theta) - A * (cos_theta - mu * sin_theta)

# Check if the block slides down. a must be positive.
# If a <= 0, the block would not slide down or the friction model is invalid.
# Our physical setup ensures a > 0.

# Step 3: Calculate the distance the block slides (L)
L = h / sin_theta

# Step 4: Calculate the time (t) using kinematics
# L = v_initial_relative * t + 0.5 * a * t^2
# Since the block starts from rest relative to the wedge, v_initial_relative = 0
# L = 0.5 * a * t^2 => t = sqrt(2 * L / a)
time_to_slide = math.sqrt(2 * L / a)

# Output the results, showing the final equation
print("This problem is solved by analyzing the forces in a non-inertial reference frame moving with the wedge.")
print(f"The calculated acceleration of the wedge is A = {A:.4f} m/s^2.")
print(f"The calculated relative acceleration of the block down the incline is a = {a:.4f} m/s^2.")
print(f"The distance the block slides along the incline is L = {L:.4f} m.")
print("\nUsing the kinematic equation L = 1/2 * a * t^2:")
print(f"The equation with the calculated values is:")
print(f"{L:.4f} m = 1/2 * {a:.4f} m/s^2 * t^2")
print(f"Solving for t = sqrt(2 * {L:.4f} / {a:.4f}) gives:")
print(f"Time for the block to reach the bottom = {time_to_slide:.4f} s.")

# Final answer in the specified format
# The output is printed above, but we also return the final numerical value.
# <<<3.5911>>>
print(f'<<<{time_to_slide}>>>')
