import math

def solve_robot_chain_problem():
    """
    Solves for the time when the chain loses contact with the ground.
    """
    # 1. Constants and derived parameters from the problem statement
    h = 1.0         # Robot height in m
    r = 0.25        # Arm length in m
    l_c = 10.0      # Chain length in m
    d = 20.0        # Visible length (diameter) of the path in m
    l_shadow = 10 * math.sqrt(3) # Shadow length of the path in m
    
    # Path radius
    R = d / 2.0
    
    # Path tilt angle (alpha)
    # l_shadow = d * cos(alpha) => cos(alpha) = l_shadow / d
    cos_alpha = l_shadow / d
    # sin(alpha) can be found using sin^2 + cos^2 = 1
    sin_alpha = math.sqrt(1 - cos_alpha**2)
    
    # 2. Derive the equation for the height of the arm tip H(t)
    # The height H(t) is composed of several parts:
    # z_center = R * sin_alpha
    # z_base_offset = R * sin(t) * sin_alpha  (where omega = 1 rad/s, so theta = t)
    # z_body_offset = h * cos_alpha
    # z_arm_offset = r * (cos(t)**2 * sin_alpha + sin(t) * cos_alpha) -> corrected derivation
    # After a detailed derivation, the height H(t) is:
    # H(t) = R*sin_alpha + R*sin(t)*sin_alpha + h*cos_alpha + r*(cos(t)**2 * sin_alpha + sin(t)*cos_alpha)
    # Let's group terms:
    # H(t) = R*sin_alpha + h*cos_alpha + (R*sin_alpha + r*cos_alpha)*sin(t) + r*sin_alpha*cos(t)**2
    # The condition for the chain losing contact is H(t) = l_c
    # This leads to a quadratic equation for s = sin(t)
    # a*s^2 + b*s + c = 0 where a=-r*sin_alpha, etc.

    # Simplified coefficients based on derivation:
    # The equation is: s^2 - (40 + sqrt(3))s + (39 - 4*sqrt(3)) = 0
    
    a = 1.0
    b = -(40 + math.sqrt(3))
    c = 39 - 4 * math.sqrt(3)

    # 3. Solve the quadratic equation for s = sin(t)
    discriminant = b**2 - 4 * a * c
    
    # The two possible solutions for s
    s1 = (-b + math.sqrt(discriminant)) / (2 * a)
    s2 = (-b - math.sqrt(discriminant)) / (2 * a)

    # 4. Find the valid solution and calculate time t
    s_valid = 0
    if -1 <= s1 <= 1:
        s_valid = s1
    elif -1 <= s2 <= 1:
        s_valid = s2
    else:
        print("Error: No physical solution found for sin(t).")
        return

    # The problem asks for the first time this happens. Since the height function
    # starts below 10m and is initially increasing, the smallest positive t is the answer.
    # This corresponds to the principal value of arcsin.
    t_solution = math.asin(s_valid)

    # 5. Output the results as requested
    print("The final equation to solve for s = sin(t) is:")
    # Printing each number in the equation: a*s^2 + b*s + c = 0
    print(f"{a} * s^2 + ({b}) * s + {c} = 0")
    print("\nSolving this equation gives two possible values for s = sin(t):")
    print(f"s1 = {s1}")
    print(f"s2 = {s2}")
    print(f"\nThe physically valid solution is s = {s_valid}")
    print("\nThe time t when the chain first loses contact with the ground is:")
    print(f"t = arcsin({s_valid}) = {t_solution} seconds")


if __name__ == "__main__":
    solve_robot_chain_problem()
