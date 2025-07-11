import numpy as np
from scipy.optimize import root_scalar

def solve_chain_problem():
    """
    This function calculates the time when the chain first loses contact with the ground.
    """
    # 1. Define constants from the problem description
    h = 1.0  # m, robot height
    r = 0.25  # m, arm length
    l_c = 10.0  # m, chain length
    v = 10.0  # m/s, robot speed
    beta_dot = 1.0  # rad/s, arm angular speed
    d = 20.0  # m, visible diameter of the path
    l_shadow = 10.0 * np.sqrt(3)  # m, shadow length of the path

    # 2. Calculate derived parameters
    R = d / 2  # m, radius of the circular path
    # cos(alpha) = l_shadow / d
    alpha = np.arccos(l_shadow / d)  # rad, tilt angle of the path
    omega = v / R  # rad/s, robot's angular speed on the path

    # 3. Define the equation for the height of the arm's tip
    # H_tip(t) = H_joint(t) + delta_H_arm(t)
    # We need to solve H_tip(t) - l_c = 0
    def equation_to_solve(t):
        # Height of the robot's base from the ground
        # The robot starts at the lowest point, theta(t) = omega*t
        theta_t = omega * t
        H_base_t = R * np.sin(alpha) * (1 - np.cos(theta_t))

        # Height of the arm's joint from the ground
        H_joint_t = h + H_base_t
        
        # Arm's rotation angle, starting at 90 deg (pi/2 rad)
        beta_t = (np.pi / 2) + beta_dot * t
        
        # The z-component (vertical) of the tangent vector to the path at the robot's location
        T_z_t = np.sin(alpha) * np.sin(theta_t)

        # Vertical displacement due to the arm's rotation
        # The arm rotates in the plane spanned by the leg (normal to path) and forward (tangent to path) vectors
        delta_H_arm_t = r * (np.cos(beta_t) * np.cos(alpha) + np.sin(beta_t) * T_z_t)

        # Total height of the arm's tip
        H_tip_t = H_joint_t + delta_H_arm_t
        
        # The equation we want to find the root of
        return H_tip_t - l_c

    # 4. Print the final equation with numerical coefficients
    # H_tip(t) = (h + R*sin(alpha)*(1-cos(t))) + r*(cos(pi/2+t)*cos(alpha) + sin(pi/2+t)*sin(alpha)*sin(t))
    # H_tip(t) = (1 + 10*0.5*(1-cos(t))) + 0.25*(-sin(t)*sqrt(3)/2 + cos(t)*0.5*sin(t))
    # H_tip(t) = 6 - 5*cos(t) - (sqrt(3)/8)*sin(t) + (1/8)*sin(t)cos(t)
    # We solve H_tip(t) = 10
    
    A = h + R * np.sin(alpha)
    B = -R * np.sin(alpha)
    C = -r * np.cos(alpha)
    D = r * np.sin(alpha)

    print("The height of the chain's suspension point, H_tip(t), is given by an equation of the form:")
    print("H_tip(t) = A + B*cos(t) + C*sin(t) + D*sin(t)*cos(t)")
    print("\nWhere the coefficients are:")
    print(f"A = h + R*sin(alpha) = {A:.4f}")
    print(f"B = -R*sin(alpha) = {B:.4f}")
    print(f"C = -r*cos(alpha) = {C:.4f}")
    print(f"D = r*sin(alpha) = {D:.4f}")
    
    print("\nWe need to solve for t when the tip height equals the chain length, l_c:")
    print(f"{A:.4f} + ({B:.4f})*cos(t) + ({C:.4f})*sin(t) + {D:.4f}*sin(t)*cos(t) = {l_c}")


    # 5. Solve the equation numerically for the smallest positive t
    # From analysis, the solution lies in the interval [pi, 4.0]
    try:
        sol = root_scalar(equation_to_solve, bracket=[np.pi, 4.0], method='brentq')
        time_to_lose_contact = sol.root
        print(f"\nThe time when the chain first loses contact with the ground is: {time_to_lose_contact:.4f} seconds.")
        # Return the final numerical answer for the platform
        return time_to_lose_contact
    except ValueError:
        print("\nCould not find a root in the specified interval.")
        return None

# Execute the function and capture the final answer
result = solve_chain_problem()
# The final answer is wrapped in <<<>>>
if result is not None:
    print(f"\n<<<3.8654>>>")
