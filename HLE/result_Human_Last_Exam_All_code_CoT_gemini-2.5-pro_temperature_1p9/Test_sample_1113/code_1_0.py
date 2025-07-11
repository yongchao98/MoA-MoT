import math

def solve_chain_problem():
    """
    Solves for the time when the chain first loses contact with the ground.
    """
    # Step 1: Define constants and calculate parameters
    h = 1.0  # m, robot height
    r = 0.25  # m, arm length
    l_c = 10.0  # m, chain length
    v = 10.0  # m/s, robot speed
    d = 20.0  # m, visible length of the path (diameter)
    l_shadow = 10 * math.sqrt(3)  # m, shadow length

    R = d / 2  # Radius of the circular path
    
    # cos(theta) = l_shadow / d
    cos_theta = l_shadow / d
    sin_theta = math.sqrt(1 - cos_theta**2)

    # Robot's angular speed on the path omega = v / R = 10 / 10 = 1 rad/s.
    # We assume phi_0=0, so phi(t) = t.
    # Arm's angular speed dot_beta = 1 rad/s. We assume beta_0=0, so beta(t) = t.
    # Let S = sin(t)

    # Step 2: Form the quadratic equation for S = sin(t)
    # The height of the chain end is given by:
    # z(t) = z_base(t) + z_head_offset + z_arm(t) - l_c
    # z_base(t) = R*sin_theta*(1 + sin(t))
    # z_head_offset = h*cos_theta
    # z_arm(t) = r*(cos^2(t)*sin_theta + sin(t)*cos_theta)
    #
    # z(t) = R*sin_theta*(1+S) + h*cos_theta + r*((1-S^2)*sin_theta + S*cos_theta) - l_c = 0
    # Grouping terms by powers of S:
    # S^2 * (-r*sin_theta) + 
    # S * (R*sin_theta + r*cos_theta) + 
    # (R*sin_theta + h*cos_theta + r*sin_theta - l_c) = 0

    A = -r * sin_theta
    B = R * sin_theta + r * cos_theta
    C = R * sin_theta + h * cos_theta + r * sin_theta - l_c

    print("The problem reduces to solving a quadratic equation for S = sin(t):")
    print(f"A*S^2 + B*S + C = 0, where:")
    print(f"A = -r * sin(theta) = {A}")
    print(f"B = R * sin(theta) + r * cos(theta) = {B}")
    print(f"C = R * sin(theta) + h * cos(theta) + r * sin(theta) - l_c = {C}\n")

    # Step 3: Solve the quadratic equation
    discriminant = B**2 - 4 * A * C

    if discriminant < 0:
        print("No real solution exists.")
        return

    S1 = (-B + math.sqrt(discriminant)) / (2 * A)
    S2 = (-B - math.sqrt(discriminant)) / (2 * A)
    
    # We are looking for the first positive time t, which means sin(t) must be positive.
    # Since sin(0) = 0 and z(0) < 0, we need the smallest positive root for S=sin(t).
    valid_solutions = []
    if -1 <= S1 <= 1:
        valid_solutions.append(S1)
    if -1 <= S2 <= 1:
        valid_solutions.append(S2)
    
    if not valid_solutions:
        print("No valid solution for sin(t) in the range [-1, 1].")
        return
        
    # We want the smallest positive S that leads to the smallest positive t
    # The chain is on the ground at t=0 (sin(t)=0).
    # The function z(S) is a downward-opening parabola, starting at z(0)=C<0.
    # It will cross z=0 at the smaller positive root.
    S = min(s for s in valid_solutions if s > 0)
    
    # Step 4: Calculate time t
    t = math.asin(S)

    print(f"Solving for S gives two roots: S1={S1}, S2={S2}")
    print(f"The valid solution for sin(t) is {S}")
    print(f"The first time t when the chain loses contact with the ground is arcsin({S}).")
    print(f"t = {t} seconds")
    
solve_chain_problem()
