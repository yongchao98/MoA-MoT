import numpy as np
from scipy.optimize import brentq

# --- Step 1: Define parameters from the problem statement ---
h = 1.0       # m, robot height
r = 0.25      # m, arm length
l_c = 10.0    # m, chain length
d = 20.0      # m, visible length of the path
l = 10 * np.sqrt(3)  # m, shadow length

# --- Step 2: Calculate geometric properties of the path ---
# Radius of the circular path
R = d / 2

# Tilt angle of the path from the horizontal
# l = d * cos(theta)
cos_theta = l / d
theta = np.arccos(cos_theta)
sin_theta = np.sin(theta)

# --- Step 3: Define the equation to solve ---
# We need to find the time 't' when the chain first loses contact with the ground.
# This occurs when the height of the arm's tip equals the length of the chain.
# The function below represents: Height_arm_tip(t) - l_c
def find_time_equation(t):
    """
    Calculates the difference between the arm tip height and the chain length.
    We are looking for the root of this function.
    """
    # Height of the robot's base above the ground
    base_height = R * sin_theta * (1 - np.cos(t))
    
    # Vertical component of the robot's leg
    leg_height_component = h * cos_theta
    
    # Vertical component of the robot's arm
    arm_height_component = r * (np.cos(t) * np.sin(t) * sin_theta + np.sin(t) * cos_theta)
    
    # Total height of the arm's tip
    total_tip_height = base_height + leg_height_component + arm_height_component
    
    return total_tip_height - l_c

# --- Step 4: Solve the equation numerically ---
# We observed that at t=pi/2 the function is negative, and at t=pi it is positive.
# This means the root lies in the interval [pi/2, pi].
# We use the Brent's method to find the root accurately.
try:
    solution_t = brentq(find_time_equation, np.pi / 2, np.pi)
    
    # --- Step 5: Print the final equation and the result ---
    print("The problem is solved by finding the root 't' for the equation: H_tip(t) = l_c")
    print("H_tip(t) = R*sin(theta)*(1 - cos(t)) + h*cos(theta) + r*sin(t)*(cos(t)*sin(theta) + cos(theta))")
    print("\nPlugging in the numerical values:")
    # The final equation to solve is:
    # 10*0.5*(1-cos(t)) + 1*(sqrt(3)/2) + 0.25*sin(t)*(cos(t)*0.5 + sqrt(3)/2) = 10
    print(f"{R} * {sin_theta:.2f} * (1 - cos(t)) + {h} * {cos_theta:.2f} + {r} * sin(t) * (cos(t) * {sin_theta:.2f} + {cos_theta:.2f}) = {l_c}")
    print("\nSolving this equation for the first time t > 0 gives the moment the chain loses contact with the ground.")
    print(f"\nThe time when the chain first loses contact with the ground is: {solution_t:.4f} seconds.")

except ValueError:
    print("Error: Could not find a solution in the specified interval.")
