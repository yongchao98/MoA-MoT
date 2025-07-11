import numpy as np
from scipy.optimize import brentq

# 1. Define all given physical parameters.
h = 1.0       # Robot height in meters
r = 0.25      # Arm length in meters (25 cm)
l_c = 10.0    # Chain length in meters
v = 10.0      # Robot speed in m/s
dot_beta = 1.0  # Arm angular speed in rad/s
d = 20.0      # Diameter of the path in meters
l_shadow = 10.0 * np.sqrt(3) # Shadow length in meters

# 2. Calculate geometric parameters of the path.
# Radius of the circular path
R = d / 2

# The tilt angle `theta` of the path is found from the shadow projection: l_shadow = d * cos(theta)
cos_theta = l_shadow / d
theta = np.arccos(cos_theta)
sin_theta = np.sin(theta)

# 3. Define the motion equations.
# Robot's angular speed `omega` = v / R
omega = v / R
# Robot's angular position alpha(t) = omega * t. Since omega=1, alpha(t) = t.
# Arm's rotation angle beta(t) = dot_beta * t. Since dot_beta=1, beta(t) = t.

# 4. Define the function for the height of the arm's tip.
# Z_tip(t) = Z_base(t) + Z_arm_component(t)
# The height of the tip of the arm is given by:
# Z_tip(t) = R*sin(theta)*(1-cos(alpha(t))) + (h + r*sin(beta(t)))*cos(theta) + r*cos(beta(t))*sin(alpha(t))*sin(theta)
def z_tip(t):
    # As omega=1 and dot_beta=1, alpha(t) = t and beta(t) = t
    alpha = t
    beta = t
    term1 = R * sin_theta * (1 - np.cos(alpha))
    term2 = (h + r * np.sin(beta)) * cos_theta
    term3 = r * np.cos(beta) * np.sin(alpha) * sin_theta
    return term1 + term2 + term3

# 5. Set up the equation to solve for the time `t`.
# The chain loses contact when Z_tip(t) = l_c.
# We need to solve the equation: z_tip(t) - l_c = 0
def equation_to_solve(t):
    return z_tip(t) - l_c

# 6. Solve the equation numerically.
# We found through analysis that the root lies between pi/2 and pi.
try:
    solution_time = brentq(equation_to_solve, np.pi/2, np.pi)
except (ValueError, RuntimeError):
    solution_time = "Error: Could not find a solution in the expected interval."

# 7. Print the final equation and the solution.
print("The equation for the height of the arm's tip (Z_tip) to equal the chain's length (l_c) is:")
print("R*sin(\u03B8)*(1-cos(t)) + (h + r*sin(t))*cos(\u03B8) + r*cos(t)*sin(t)*sin(\u03B8) = l_c")
print("\nSubstituting the numerical values:")

final_equation_str = (
    f"{R} * {sin_theta:.3f} * (1 - cos(t)) + "
    f"({h} + {r} * sin(t)) * {cos_theta:.3f} + "
    f"{r} * cos(t) * sin(t) * {sin_theta:.3f} = {l_c}"
)
print(final_equation_str)

print(f"\nSolving this equation for t gives the time when the chain first loses contact with the ground.")
print(f"Time (t) = {solution_time:.4f} seconds.")
