import math

# Define the given parameters
m = 0.1  # mass of the block in kg
M = 10.0 # mass of the wedge in kg
h = 2.0  # height of the wedge in m
theta_deg = 30.0 # angle of the wedge in degrees
mu = 0.5 # coefficient of kinetic friction
g = 10.0 # acceleration due to gravity in m/s^2

# Convert the angle from degrees to radians for math functions
theta = math.radians(theta_deg)

# Pre-calculate sin and cos of the angle
sin_theta = math.sin(theta)
cos_theta = math.cos(theta)

# The equations for the system are derived from Newton's second law.
# 1) Perpendicular to incline (on block): N - mg*cos(th) + m*A_x*sin(th) = 0
# 2) Along incline (on block): mg*sin(th) - mu*N - m*A_x*cos(th) = m*a_rel
# 3) Horizontal (on wedge): N*sin(th) + mu*N*cos(th) = M*A_x

# Solving this system of equations for A_x (acceleration of the wedge) and a_rel (relative acceleration of the block)

# Calculate the horizontal acceleration of the wedge, A_x
ax_numerator = m * g * cos_theta * (sin_theta + mu * cos_theta)
ax_denominator = M + m * (sin_theta**2) + m * mu * sin_theta * cos_theta
A_x = ax_numerator / ax_denominator

# Calculate the acceleration of the block relative to the wedge, a_rel
term1 = g * (sin_theta - mu * cos_theta)
term2 = A_x * (cos_theta - mu * sin_theta)
a_rel = term1 - term2

# Calculate the distance 'L' the block must slide down the incline
L = h / sin_theta

# Calculate the time 't' using the kinematic equation: L = 0.5 * a_rel * t^2
# The block starts from rest relative to the wedge, so initial velocity is 0.
time_squared = 2 * L / a_rel
t = math.sqrt(time_squared)

# Print the final equation and the values of its components as requested
print("The final equation for time 't' is derived from the kinematic formula L = (1/2) * a_rel * t^2")
print("Solving for t gives: t = sqrt(2 * L / a_rel)")
print(f"The distance the block slides, L = {L:.4f} m")
print(f"The block's acceleration relative to the wedge, a_rel = {a_rel:.4f} m/s^2")
print(f"The time it takes for the block to reach the bottom is t = sqrt(2 * {L:.4f} / {a_rel:.4f}) = {t:.4f} s")

# Final answer in the specified format
print(f"\n<<<{t}>>>")