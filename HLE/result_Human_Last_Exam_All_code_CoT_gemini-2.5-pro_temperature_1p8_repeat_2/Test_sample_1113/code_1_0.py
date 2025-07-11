import numpy as np
from scipy.optimize import fsolve

# Step 1: Define given constants
h = 1.0       # Robot height in m
r = 0.25      # Arm length in m
l_c = 10.0    # Chain length in m
d = 20.0      # Visible diameter of the path in m
l_shadow = 10 * np.sqrt(3) # Shadow length of the path in m
v = 10.0      # Robot speed in m/s
dot_beta = 1.0 # Arm angular speed in rad/s

# Step 2: Calculate derived parameters
R = d / 2.0  # Radius of the circular path
omega = v / R # Angular speed of the robot on the path

# Calculate the tilt angle theta of the path
cos_theta = l_shadow / d
sin_theta = np.sqrt(1 - cos_theta**2)

# Step 3: Define the function for the height of the arm tip z_tip(t)
# We need to solve z_tip(t) - l_c = 0.
# The assumption is that the robot starts at alpha_0 = 0.
# So alpha(t) = omega*t. Since omega=1, alpha(t) = t.
# The arm rotates with dot_beta = 1, and starts 'forward', perpendicular to the leg.
# This means beta(t) = np.pi/2 + t
def get_z_tip(t, R, h, r, cos_theta, sin_theta):
    """Calculates the vertical height of the arm tip at time t."""
    # alpha(t) = t, beta(t) = np.pi/2 + t
    alpha = t
    
    # Height of the robot's feet on the circular path
    z_base = R * sin_theta * (1 + np.sin(alpha))
    
    # Vertical height from the robot's body
    z_body = h * cos_theta
    
    # Vertical component from the rotating arm
    # The arm vector in the tilted frame (x',y',z') has components:
    # v_arm_y_prime = -r * np.sin(beta) * np.cos(alpha)
    # v_arm_z_prime = r * np.cos(beta)
    # The transformation to the ground frame z-component is:
    # z_arm = v_arm_y_prime * sin_theta + v_arm_z_prime * cos_theta
    beta = np.pi/2 + t
    z_arm = -r * np.sin(beta) * np.cos(alpha) * sin_theta + r * np.cos(beta) * cos_theta
    
    return z_base + z_body + z_arm

# The equation to solve is f(t) = 0
def equation_to_solve(t):
    return get_z_tip(t, R, h, r, cos_theta, sin_theta) - l_c

# Step 4: Print the equation and solve it
print("Solving for the time 't' when the chain loses contact with the ground.")
print("This occurs when the height of the arm tip, z_tip(t), equals the chain length, l_c.")
print("\n--- Given and Derived Parameters ---")
print(f"Robot height (h) = {h} m")
print(f"Arm length (r) = {r} m")
print(f"Chain length (l_c) = {l_c} m")
print(f"Path radius (R) = {R} m")
print(f"Path tilt angle (theta) = {np.rad2deg(np.arccos(cos_theta)):.2f} degrees")

print("\n--- The Equation ---")
print("z_tip(t) = l_c")
print("R*sin(theta)*(1+sin(t)) + h*cos(theta) - r*cos(t)*sin(t+pi/2)*sin(theta) + r*cos(t+pi/2)*cos(theta) = l_c")
print("Substituting numerical values, we solve the following equation for t:")

# Equation with numerical values for printing
R_sin_theta = R * sin_theta
h_cos_theta = h * cos_theta
r_cos_theta = r * cos_theta
r_sin_theta = r * sin_theta

# Simplified equation printout
# 10 = R*sin(th) + h*cos(th) + R*sin(th)*sin(t) - r*cos(th)*sin(t) - r*sin(th)*cos^2(t)
final_eq_str = (f"{l_c:.4f} = {R_sin_theta:.4f}*(1 + sin(t)) + {h_cos_theta:.4f} + "
                f"{r:.4f}*(-sin(t)*{cos_theta:.4f} - cos(t)^2*{sin_theta:.4f})")

# A more readable rearrangement for printing:
# A + B*sin(t) - C*cos^2(t) - l_c = 0
A_term = R_sin_theta + h_cos_theta - l_c
B_term = R_sin_theta - r * cos_theta
C_term = r * sin_theta

print(f"{A_term:.4f} + ({B_term:.4f})*sin(t) - ({C_term:.4f})*cos(t)^2 = 0")


# Step 5: Solve the equation numerically
# From plotting or analyzing the function, the first positive root is near t=1.0
initial_guess = 1.0
solution_t = fsolve(equation_to_solve, initial_guess)

print("\n--- Solution ---")
# fsolve returns an array, so we take the first element
time_to_lose_contact = solution_t[0]
print(f"The chain first loses contact with the ground at t = {time_to_lose_contact:.4f} seconds.")

# Final answer format
print(f"\nFinal numerical answer in the required format:")
print(f"<<<{time_to_lose_contact:.4f}>>>")
