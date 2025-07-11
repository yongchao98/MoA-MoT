import math

def solve_chain_contact_problem():
    """
    Solves for the time when the chain first loses contact with the ground.
    """
    # Given parameters from the problem description
    h = 1.0       # Robot height in meters
    r = 0.25      # Arm length in meters (25 cm)
    l_c = 10.0    # Chain length in meters
    v = 10.0      # Robot speed in m/s
    d = 20.0      # Visible diameter of the path in meters
    l_shadow = 10 * math.sqrt(3) # Shadow length in meters

    # Step 1: Determine the geometry of the circular path
    R = d / 2.0
    # From l_shadow = d * cos(alpha), we find the tilt angle alpha
    cos_alpha = l_shadow / d
    alpha = math.acos(cos_alpha)
    sin_alpha = math.sin(alpha)

    # Step 2: Formulate the height equation for the end of the chain
    # The chain loses contact when its end is at ground level, i.e., z_end(t) = 0.
    # The height z_end(t) can be expressed as a function of s = sin(t).
    # z_end(s) = a*s^2 + b*s + c = 0
    # where the coefficients a, b, and c are derived from the system's geometry and kinematics.
    
    # Derivation of coefficients:
    # z_end(t) = z_robot(t) + z_leg(t) + z_arm(t) - l_c
    # z_robot(t) = R * sin(alpha) * (1 + sin(t))
    # z_leg(t) = h * cos(alpha)
    # z_arm(t) = r * (sin(alpha)*cos^2(t) + cos(alpha)*sin(t))
    # Substituting s = sin(t) and cos^2(t) = 1 - s^2, we get:
    # z_end(s) = R*sin(a)*(1+s) + h*cos(a) + r*(sin(a)*(1-s^2) + cos(a)*s) - l_c
    # Grouping terms by powers of s:
    # z_end(s) = (-r*sin(a))*s^2 + (R*sin(a) + r*cos(a))*s + (R*sin(a) + h*cos(a) + r*sin(a) - l_c)
    
    a_quad = -r * sin_alpha
    b_quad = R * sin_alpha + r * cos_alpha
    c_quad = (R * sin_alpha) + (h * cos_alpha) + (r * sin_alpha) - l_c

    print("The height of the chain's end is given by z(t). The chain loses contact when z(t) = 0.")
    print("This condition forms a quadratic equation in terms of s = sin(t):")
    print(f"a*s^2 + b*s + c = 0")
    print("The coefficients are:")
    print(f"a = -r * sin(alpha) = {a_quad}")
    print(f"b = R * sin(alpha) + r * cos(alpha) = {b_quad}")
    print(f"c = R*sin(alpha) + h*cos(alpha) + r*sin(alpha) - l_c = {c_quad}")
    print("-" * 30)

    # Step 3: Solve the quadratic equation for s = sin(t)
    discriminant = b_quad**2 - 4 * a_quad * c_quad

    if discriminant < 0:
        print("No real solution exists. The chain never reaches ground level.")
        return

    s1 = (-b_quad + math.sqrt(discriminant)) / (2 * a_quad)
    s2 = (-b_quad - math.sqrt(discriminant)) / (2 * a_quad)

    print(f"Solving the quadratic equation gives two possible values for s:")
    print(f"s1 = {s1}")
    print(f"s2 = {s2}")
    
    # Step 4: Find the smallest positive time t
    # A valid solution for s must be in the range [-1, 1]
    valid_s = []
    if -1 <= s1 <= 1:
        valid_s.append(s1)
    if -1 <= s2 <= 1:
        valid_s.append(s2)

    if not valid_s:
        print("No valid solution for s = sin(t) in the range [-1, 1].")
        return

    # For each valid s, find the smallest positive t such that sin(t) = s
    possible_t = []
    for s in valid_s:
        # The principal value from asin is in [-pi/2, pi/2]
        t_asin = math.asin(s)
        # If t_asin is positive, it's the smallest positive solution.
        if t_asin > 0:
            possible_t.append(t_asin)
        # Other potential positive solutions are pi - t_asin, 2*pi + t_asin, etc.
        # We only need the first time this happens.
        
    final_t = min(possible_t)

    print("-" * 30)
    print("The physically valid solution for s = sin(t) is the one between -1 and 1.")
    print(f"The first time the chain loses contact corresponds to the smallest positive t.")
    print(f"The time is t = arcsin({min(valid_s):.4f})")
    print(f"t = {final_t:.4f} seconds")
    return final_t

# Run the solver and store the final answer
result = solve_chain_contact_problem()
print(f"\n<<< {result:.4f} >>>")
