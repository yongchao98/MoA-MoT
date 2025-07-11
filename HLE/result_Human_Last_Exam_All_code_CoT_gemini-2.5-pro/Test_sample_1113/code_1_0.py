import numpy as np

def solve_robot_chain_problem():
    """
    This function calculates the time when the robot's chain first loses contact with the ground.
    """
    # Define problem parameters from the user's request
    h = 1.0       # robot height in meters
    r = 0.25      # arm length in meters
    l_c = 10.0    # chain length in meters
    v = 10.0      # robot speed in m/s
    beta_dot = 1.0  # arm angular speed in rad/s
    d = 20.0      # path diameter in meters
    l_shadow = 10 * np.sqrt(3) # path shadow length in meters

    # --- Step 1: Determine the geometry of the tilted path ---
    # The diameter of the circular path is d.
    R = d / 2.0

    # The path is tilted by an angle alpha. The sun is overhead.
    # The shadow is an ellipse with major axis d and minor axis d*cos(alpha).
    # The problem states the shadow length is l_shadow, which corresponds to the minor axis.
    # So, d * cos(alpha) = l_shadow
    cos_alpha = l_shadow / d
    sin_alpha = np.sqrt(1 - cos_alpha**2)

    # --- Step 2: Formulate the height equation for the chain's tip ---
    # The robot's angular position on the circle is theta(t) = (v/R)*t = (10/10)*t = t.
    # The arm's rotational angle is also a function of time, gamma(t) = beta_dot*t = t,
    # based on the initial conditions.
    # The chain loses contact when its tip's height Z_tip(t) = 0 for the first time (for t > 0).
    # The full equation for the height is:
    # Z_tip(t) = R*sin_alpha*(1+sin(t)) + h*cos_alpha + r*(cos^2(t)*sin_alpha - sin(t)*cos_alpha) - l_c = 0
    # This can be rewritten as a quadratic equation for x = sin(t): a*x^2 + b*x + c = 0

    # --- Step 3: Define coefficients and solve the equation ---
    a_coeff = -r * sin_alpha
    b_coeff = (R * sin_alpha) - (r * cos_alpha)
    c_coeff = (R * sin_alpha) + (h * cos_alpha) + (r * sin_alpha) - l_c

    # Solve the quadratic equation using the quadratic formula
    discriminant = b_coeff**2 - 4 * a_coeff * c_coeff
    
    if discriminant < 0:
        print("Error: No real solution exists for the time.")
        return

    # Calculate the two possible solutions for x = sin(t)
    x1 = (-b_coeff + np.sqrt(discriminant)) / (2 * a_coeff)
    x2 = (-b_coeff - np.sqrt(discriminant)) / (2 * a_coeff)

    # The physically valid solution for sin(t) must be in the range [-1, 1]
    sin_t_sol = None
    if -1 <= x1 <= 1:
        sin_t_sol = x1
    elif -1 <= x2 <= 1:
        sin_t_sol = x2

    if sin_t_sol is None:
        print("Error: No physically valid solution for sin(t) was found.")
        return

    # The chain starts on the ground, so we seek the smallest positive time t.
    # This corresponds to the principal value of arcsin.
    time = np.arcsin(sin_t_sol)

    print("The final equation for the height of the chain tip Z(t) is set to 0:")
    print("Z(t) = R*sin(alpha)*(1+sin(t)) + h*cos(alpha) + r*(cos^2(t)*sin(alpha) - sin(t)*cos(alpha)) - l_c = 0")
    print("\nSubstituting the numerical values, the equation with each term calculated is:")
    
    # Output each number/term in the final equation
    term1_val = R * sin_alpha
    term2_val = h * cos_alpha
    term3_val_cos2 = r * sin_alpha
    term3_val_sin = r * cos_alpha
    print(f"{term1_val:.4f} * (1 + sin(t)) + {term2_val:.4f} + {term3_val_cos2:.4f} * cos^2(t) - {term3_val_sin:.4f} * sin(t) - {l_c:.4f} = 0")

    print("\nThis simplifies to the quadratic equation a*x^2 + b*x + c = 0, where x = sin(t):")
    print(f"{a_coeff:.4f}*x^2 + {b_coeff:.4f}*x + {c_coeff:.4f} = 0")
    
    print(f"\nThe physically valid solution is sin(t) = {sin_t_sol:.4f}")
    
    print(f"\nThe time 't' when the chain first loses contact with the ground is:")
    print(f"t = arcsin({sin_t_sol:.4f}) = {time:.4f} seconds")

solve_robot_chain_problem()