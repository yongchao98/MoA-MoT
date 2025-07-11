import math

def solve_sliding_time():
    """
    Calculates the time it takes for a block to slide down a moving wedge with friction.
    """
    # Given parameters
    m = 0.1  # kg, mass of the block
    M = 10.0  # kg, mass of the wedge
    h = 2.0  # m, initial height of the block
    theta_deg = 30.0  # degrees, angle of the wedge
    mu = 0.5  # coefficient of kinetic friction
    g = 10.0  # m/s^2, acceleration due to gravity
    # v0 = 1.0 m/s is the initial velocity of the system. It does not affect the
    # relative acceleration or the time to slide down, as the block starts
    # from rest relative to the wedge.

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # Step 1: Calculate the distance the block slides along the wedge's surface
    s = h / sin_theta
    print(f"The angle of the wedge is {theta_deg} degrees.")
    print(f"The block starts at a height h = {h} m.")
    print(f"The distance the block needs to slide is s = h / sin({theta_deg}) = {h:.4f} / {sin_theta:.4f} = {s:.4f} m.")
    print("-" * 30)

    # Step 2: Calculate the acceleration of the block relative to the wedge (a_rel)
    # The formula is derived from Newton's second law in a non-inertial reference frame.
    # a_rel = g * ((M-m)*sin(theta) - mu*(M+m)*cos(theta)) / (M - m*sin(theta)**2 - mu*m*sin(theta)*cos(theta))
    
    a_rel_numerator = g * ((M - m) * sin_theta - mu * (M + m) * cos_theta)
    a_rel_denominator = M - m * sin_theta**2 - mu * m * sin_theta * cos_theta
    a_rel = a_rel_numerator / a_rel_denominator
    
    print("The acceleration of the block relative to the wedge (a_rel) is calculated.")
    print(f"Numerator of a_rel = g*((M-m)*sin(theta) - mu*(M+m)*cos(theta)) = {a_rel_numerator:.4f} m^2/s^2")
    print(f"Denominator of a_rel = M - m*sin(theta)^2 - mu*m*sin(theta)*cos(theta) = {a_rel_denominator:.4f} kg")
    print(f"Relative acceleration a_rel = {a_rel_numerator:.4f} / {a_rel_denominator:.4f} = {a_rel:.4f} m/s^2.")
    print("-" * 30)

    # Step 3: Calculate the time it takes for the block to slide down
    # Using the kinematic equation: s = (1/2) * a_rel * t^2
    # Solving for t: t = sqrt(2 * s / a_rel)
    if a_rel <= 0:
        print("The relative acceleration is zero or negative. The block will not slide down.")
        time = float('inf')
    else:
        time = math.sqrt(2 * s / a_rel)

    print("The time to slide down is calculated using the equation: t = sqrt(2 * s / a_rel)")
    print("Plugging in the numbers:")
    # The requirement is to output each number in the final equation.
    print(f"t = sqrt(2 * {s:.4f} / {a_rel:.4f})")
    print(f"t = sqrt({2 * s:.4f} / {a_rel:.4f})")
    print(f"t = sqrt({(2 * s / a_rel):.4f})")
    print(f"t = {time:.4f} s")
    
    return time

# Execute the function and store the final answer.
final_time = solve_sliding_time()
print(f"<<<{final_time:.4f}>>>")
