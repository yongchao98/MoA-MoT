import numpy as np

# 1. Define constants
h = 1.0  # m, robot height
r = 0.25  # m, arm length
lc = 10.0  # m, chain length
v = 10.0  # m/s, robot speed
d = 20.0  # m, visible path length (diameter)
l = 10 * np.sqrt(3)  # m, shadow length
dot_beta = 1.0 # rad/s, arm angular speed

# 2. Deconstruct Geometry
R = d / 2.0  # m, radius of the circular path
cos_theta = l / d
theta = np.arccos(cos_theta)
sin_theta = np.sin(theta)

# 3. Define robot's and arm's motion parameters
omega = v / R  # Robot's angular speed on the path, rad/s
# From the problem statement analysis, we assume the starting angular position
# alpha_start = pi/2, so alpha(t) = pi/2 + omega*t
# Arm angle beta(t) = dot_beta*t
# Note: omega = 10/10 = 1 and dot_beta = 1, so alpha(t) = pi/2 + t and beta(t) = t.
# The variable 't' represents both the time and the angle for beta and the change in alpha.

# 4. Formulate and solve the equation Z_tip(t) = lc
# The height equation Z_tip(t) is:
# R*sin(theta)*(1 - cos(alpha(t))) + h*cos(theta) + r*(sin(alpha(t))*cos(beta(t))*sin(theta) + sin(beta(t))*cos(theta)) = lc
# Substituting alpha(t) = pi/2 + t and beta(t) = t:
# cos(pi/2 + t) = -sin(t)
# sin(pi/2 + t) = cos(t)
# The equation becomes:
# R*sin(theta)*(1 + sin(t)) + h*cos(theta) + r*(cos(t)*cos(t)*sin(theta) + sin(t)*cos(theta)) = lc
# R*sin(theta) + R*sin(theta)*sin(t) + h*cos(theta) + r*sin(theta)*cos^2(t) + r*cos(theta)*sin(t) = lc
# Substitute cos^2(t) = 1 - sin^2(t):
# R*sin(theta) + R*sin(theta)*sin(t) + h*cos(theta) + r*sin(theta)*(1-sin^2(t)) + r*cos(theta)*sin(t) = lc
# This can be rearranged into a quadratic equation for x = sin(t): Ax^2 + Bx + C = 0

# Coefficients of the quadratic equation A*x^2 + B*x + C = 0 for x = sin(t)
A = -r * sin_theta
B = R * sin_theta + r * cos_theta
C = R * sin_theta + h * cos_theta + r * sin_theta - lc

# 5. Plug in the values and calculate the final answer
# Calculate coefficients with the given numbers
A_val = -r * sin_theta
B_val = R * sin_theta + r * cos_theta
C_val = R * sin_theta + h * cos_theta + r * sin_theta - lc

# Solve the quadratic equation: A*x^2 + B*x + C = 0
roots = np.roots([A_val, B_val, C_val])

# The valid root for sin(t) must be between -1 and 1
valid_sin_t = 0
for root in roots:
    if -1 <= root <= 1 and np.isreal(root):
        valid_sin_t = np.real(root)
        break

# Calculate the time t
# We want the smallest positive time, so we take the principal arcsin value
time = np.arcsin(valid_sin_t)

print(f"Radius of circular path (R): {R} m")
print(f"Tilt angle (theta): {np.degrees(theta):.2f} degrees")
print(f"Equation for sin(t): {A_val:.4f}*sin(t)^2 + {B_val:.4f}*sin(t) + {C_val:.4f} = 0")
print(f"Solving for sin(t), we get: sin(t) = {valid_sin_t:.4f}")
print(f"The chain will first lose contact with the ground at t = {time:.4f} seconds.")
<<<0.9008>>>