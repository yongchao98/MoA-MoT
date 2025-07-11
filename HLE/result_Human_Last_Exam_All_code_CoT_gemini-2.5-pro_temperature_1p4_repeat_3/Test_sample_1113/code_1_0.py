import numpy as np

def solve_robot_chain_problem():
    """
    Calculates the time when the chain first loses contact with the ground.
    """
    # Step 1: Define all constants from the problem description
    h = 1.0       # m, robot height
    r = 0.25      # m, arm length (25 cm)
    lc = 10.0     # m, chain length
    v = 10.0      # m/s, robot speed
    d = 20.0      # m, diameter of the path
    l_shadow = 10.0 * np.sqrt(3) # m, shadow length
    dot_beta = 1.0 # rad/s, arm angular speed

    # Step 2: Calculate geometric and motion parameters
    R = d / 2  # m, radius of the circular path
    
    # The tilt angle alpha is found from the shadow projection: l_shadow = d * cos(alpha)
    cos_alpha = l_shadow / d
    # Using the identity sin^2(alpha) + cos^2(alpha) = 1
    sin_alpha = np.sqrt(1 - cos_alpha**2)
    
    # Robot's angular speed on the path
    omega = v / R 

    # Since omega = 1 rad/s, the robot's angle on the path is simply t.

    # Step 3: Formulate the equation for the z-coordinate of the chain's end
    # We need to find t such that z_chain_end(t) = 0.
    # z_chain_end(t) = z_tip(t) - lc = 0  =>  z_tip(t) = lc
    
    # The z-coordinate of the tip of the arm is composed of several parts:
    # 1. Height of the circle's center: h_center = R * sin_alpha
    # 2. Vertical component from robot's position on the circle: R * sin(t) * sin_alpha
    # 3. Vertical component from robot's own height: h * cos_alpha
    # 4. Vertical component from the rotating arm: r * (cos^2(t)*sin_alpha + sin(t)*cos_alpha)
    
    # Let's set up the equation z_tip(t) = lc:
    # (R*sin_alpha) + (R*sin(t)*sin_alpha) + (h*cos_alpha) + r*(np.cos(t)**2*sin_alpha + np.sin(t)*cos_alpha) = lc
    
    # Substitute cos^2(t) = 1 - sin^2(t) and let s = sin(t).
    # The equation becomes a quadratic equation of the form a*s^2 + b*s + c = 0.
    
    h_center = R * sin_alpha
    
    # Coefficient for s^2
    a = -r * sin_alpha
    
    # Coefficient for s
    b = R * sin_alpha + r * cos_alpha
    
    # Constant term
    c = h_center + h * cos_alpha + r * sin_alpha - lc

    # We now have the equation: a*s^2 + b*s + c = 0
    # For clarity in the output, we multiply by -1/a to make the s^2 coefficient 1.
    B_eq = b / -a
    C_eq = c / -a
    
    # Step 4: Solve the quadratic equation
    discriminant = b**2 - 4 * a * c
    
    # There are two potential solutions for s
    s1 = (-b + np.sqrt(discriminant)) / (2 * a)
    s2 = (-b - np.sqrt(discriminant)) / (2 * a)
    
    # We choose the physically valid solution for s = sin(t), which must be in the range [-1, 1].
    # At t=0, the chain starts on the ground, so z_chain_end(0) <= 0. As t increases, z will increase.
    # The first time it leaves the ground corresponds to the smallest positive t.
    if -1 <= s1 <= 1:
        s = s1
    else:
        s = s2
        
    # Step 5: Calculate the time t
    # The first time (smallest positive t) is given by arcsin(s)
    t = np.arcsin(s)
    
    # Print the explanation and the final result
    print("To find when the chain first loses contact with the ground, we solve for the time 't' where the height of the chain's end is zero.")
    print("This leads to a quadratic equation for s = sin(t):")
    print("\n  A * sin(t)^2 + B * sin(t) + C = 0")
    print("\nAfter substituting the given values and simplifying, the equation becomes:")
    print(f"  (1.000) * sin(t)^2 + ({B_eq:.3f}) * sin(t) + ({C_eq:.3f}) = 0")
    print("\nSolving this equation for sin(t) gives a valid root:")
    print(f"  sin(t) = {s:.4f}")
    print("\nThe first time 't' (for t > 0) that satisfies this is:")
    print(f"  t = arcsin({s:.4f})")

solve_robot_chain_problem()

# Final Answer Calculation
h = 1.0; r = 0.25; lc = 10.0; v = 10.0; d = 20.0; l_shadow = 10.0 * np.sqrt(3)
R = d / 2
cos_alpha = l_shadow / d
sin_alpha = np.sqrt(1 - cos_alpha**2)
h_center = R * sin_alpha
a = -r * sin_alpha
b = R * sin_alpha + r * cos_alpha
c = h_center + h * cos_alpha + r * sin_alpha - lc
discriminant = b**2 - 4 * a * c
s1 = (-b + np.sqrt(discriminant)) / (2 * a)
s2 = (-b - np.sqrt(discriminant)) / (2 * a)
if -1 <= s1 <= 1: s = s1
else: s = s2
t = np.arcsin(s)
print(f"\n  t = {t:.4f} seconds")