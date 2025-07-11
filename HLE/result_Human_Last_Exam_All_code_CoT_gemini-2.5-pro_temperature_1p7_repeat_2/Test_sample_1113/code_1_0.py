import numpy as np
from scipy.optimize import root_scalar

# Step 1: Define constants from the problem description
h = 1.0        # robot height in meters
r = 0.25       # arm length in meters
l_c = 10.0     # chain length in meters
v = 10.0       # robot speed in m/s
beta_dot = 1.0 # arm angular speed in rad/s
d = 20.0       # visible path length (diameter) in meters
l_shadow = 10 * np.sqrt(3) # shadow length in meters

# Step 2: Calculate geometric and motion parameters
R = d / 2      # Radius of the circular path
omega = v / R  # Angular velocity of the robot on the path

# Calculate the tilt angle alpha
# l_shadow = d * cos(alpha)
cos_alpha = l_shadow / d
alpha = np.arccos(cos_alpha)

# The center of the circle is at height Z_center = R * sin(alpha)
Z_center = R * np.sin(alpha)

# We assume the robot starts at the lowest point of the circle (theta_0 = pi).
# The position angle is theta_pos(t) = omega*t + pi = t + pi.
# The arm rotation angle is beta(t) = beta_dot*t = t.

# Step 3: Define the function for the height of the arm tip
# The height is H_tip(t) = Z_center - x_local_tip(t)*sin(alpha) + z_local_tip(t)*cos(alpha)
# where x_local_tip(t) = R*cos(t+pi) + r*cos(t)*(-sin(t+pi))
# and z_local_tip(t) = h + r*sin(t)
# After simplification, we get the equation H_tip(t) - l_c = 0.

def height_equation(t):
    """
    Represents the equation H_tip(t) - l_c = 0.
    The function returns H_tip(t) - l_c.
    """
    # H_tip(t) = Z_center + R*np.sin(alpha)*np.cos(t) - (r/2)*np.sin(t)*np.cos(t)*np.sin(alpha) \
    #            + h*np.cos(alpha) + r*np.sin(t)*np.cos(alpha)
    # This simplifies to the expression below.
    
    term1 = 5 * np.cos(t)
    term2 = - (r / 2) * np.sin(2 * t) * np.sin(alpha) # Using sin(t)cos(t) = sin(2t)/2
    term3 = h * np.cos(alpha)
    term4 = r * np.sin(t) * np.cos(alpha)
    
    # H_tip(t) = Z_center + term1 - ... This is also complex. Let's use the final clear form.
    # 5*cos(t) - 0.125*sin(t)*cos(t) + (sqrt(3)/8)*sin(t) = 5 - sqrt(3)/2
    
    lhs = 5 * np.cos(t) - 0.125 * np.sin(t) * np.cos(t) + (np.sqrt(3) / 8) * np.sin(t)
    rhs = 5 - np.sqrt(3) / 2
    
    return lhs - rhs

# Step 4: Solve the equation and print the output
# Based on analysis, a root exists between 0 and pi/2.
# We will use a numerical solver to find it.

# Coefficients for the equation A*cos(t) + B*sin(t)*cos(t) + C*sin(t) = D
A = 5.0
B = -0.125
C = np.sqrt(3) / 8
D = 5 - np.sqrt(3) / 2

print("The problem reduces to finding the smallest positive time 't' that solves the following equation for the height of the chain's attachment point being equal to the chain's length:")
print("H_tip(t) = l_c")
print("\nThis results in the transcendental equation:")
print(f"({A:.4f})*cos(t) + ({B:.4f})*sin(t)*cos(t) + ({C:.4f})*sin(t) = {D:.4f}")
print("\nSolving this equation numerically...")

# Find the smallest positive root
try:
    sol = root_scalar(height_equation, bracket=[0, np.pi/2], method='brentq')
    if sol.converged:
        time_to_lose_contact = sol.root
        print(f"\nThe chain will first lose contact with the ground at t = {time_to_lose_contact:.4f} seconds.")
    else:
        print("\nNumerical solver did not converge to a solution.")
except ValueError:
    print("\nCould not find a root in the specified interval. The function may not cross zero.")

# Final answer to be encapsulated
final_answer = root_scalar(height_equation, bracket=[0, np.pi / 2], method='brentq').root
print(f"\n<<<%.4f>>>" % final_answer)
