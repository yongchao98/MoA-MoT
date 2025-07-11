import numpy as np
from scipy.optimize import fsolve

# Define the constants from the problem description
h = 1.0  # m, robot height
r = 0.25  # m, arm length
l_c = 10.0  # m, chain length
d = 20.0  # m, path diameter
R = d / 2.0  # m, path radius
shadow_l = 10 * np.sqrt(3) # m, shadow length
v = 10.0  # m/s, robot speed
beta_dot = 1.0  # rad/s, arm angular speed

# Calculate geometric and motion parameters
# The tilt angle of the path
theta = np.arccos(shadow_l / d)
# Angular speed of the robot on the path
omega_path = v / R

# The problem states the robot reaches the highest point during the first quarter of its journey.
# This implies the starting angle phi_0 is in the range (3*pi/2, 2*pi).
# In such problems, it's common to assume a boundary condition for a unique solution.
# We'll assume the robot starts at the lowest possible point in this range.
phi_0 = 3 * np.pi / 2

# The height of the chain end must be zero for it to lose contact with the ground.
# This means the height of the arm's tip must equal the chain's length.
# Z_tip(t) = l_c
# The equation for Z_tip is derived from the geometry:
# Z_tip = Z_center + Z_path_component + Z_robot_body_component + Z_arm_component
# Z_tip(t) = (R*np.sin(theta)) + (R*np.cos(phi(t))*np.sin(theta)) + (h*np.cos(theta)) + Z_arm(t)
# where phi(t) = phi_0 + omega_path * t and Z_arm is the vertical contribution from the arm's orientation.
# Z_arm(t) = r * (np.cos(beta(t))*np.cos(theta) - np.sin(beta(t))*np.sin(phi(t))*np.sin(theta))
# where beta(t) = pi/2 + beta_dot * t
# Setting Z_tip(t) = l_c gives us the equation to solve for t.

def equation_to_solve(t):
    """
    This function represents the condition Z_tip(t) - l_c = 0.
    We are looking for the root of this function.
    """
    # Robot's angle on the path at time t
    phi = phi_0 + omega_path * t
    # Arm's angle at time t
    beta = np.pi / 2 + beta_dot * t
    
    # Z-coordinate of the center of the circular path
    z_center = R * np.sin(theta)
    
    # Height of the shoulder joint
    z_shoulder = z_center + R * np.cos(phi) * np.sin(theta) + h * np.cos(theta)
    
    # Vertical component of the arm vector
    z_arm = r * (np.cos(beta) * np.cos(theta) - np.sin(beta) * np.sin(phi) * np.sin(theta))
    
    # Height of the arm's tip
    z_tip = z_shoulder + z_arm
    
    # The equation to find the root of
    return z_tip - l_c

# Use a numerical solver to find the smallest positive time 't'
# We start the search from a small positive value to find the first time it loses contact.
solution_time = fsolve(equation_to_solve, x0=0.1)[0]

# Print the final result in a descriptive way
print("The robot's parameters are:")
print(f"Height h = {h} m")
print(f"Arm length r = {r} m")
print(f"Chain length l_c = {l_c} m")
print(f"Path radius R = {R} m")
print(f"Path tilt angle theta = {np.degrees(theta):.2f} degrees")
print(f"Robot linear speed v = {v} m/s")
print(f"Arm angular speed beta_dot = {beta_dot} rad/s")
print("\nThe condition for the chain to lose contact with the ground is when the height of the arm's tip equals the length of the chain.")
print(f"Z_tip(t) = {l_c:.2f}")
print("\nSolving the equation Z_tip(t) - l_c = 0 for the smallest positive time t:")
print(f"The chain will first lose contact with the ground at approximately t = {solution_time:.3f} seconds.")
