import math

def solve():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    """
    # Given parameters
    m = 0.1  # kg, mass of the block
    M = 10.0 # kg, mass of the wedge
    h = 2.0  # m, height of the wedge
    theta_deg = 30.0 # degrees, angle of the wedge
    mu = 0.5 # coefficient of kinetic friction
    g = 10.0 # m/s^2, acceleration due to gravity
    # v0 = 1.0 m/s is the initial velocity of the system. It does not affect
    # the forces or accelerations, and thus does not affect the time to slide.

    # Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # --- Calculation from derived equations of motion ---

    # We solve a system of equations derived from Newton's second law for the
    # block and the wedge. This allows us to find the acceleration of the wedge (A_x)
    # and the acceleration of the block relative to the wedge (a_rel).
    #
    # The derived expression for A_x is:
    # A_x = (m*g*(mu*cos_theta**2 - sin_theta*cos_theta)) / (M + m*(mu*sin_theta*cos_theta - sin_theta**2))
    
    # Let's calculate the terms for A_x
    num_Ax = m * g * (mu * cos_theta**2 - sin_theta * cos_theta)
    den_Ax = M + m * (mu * sin_theta * cos_theta - sin_theta**2)
    
    if den_Ax == 0:
        A_x = 0.0 # Avoid division by zero, not expected in this problem
    else:
        A_x = num_Ax / den_Ax

    # The derived expression for a_rel is:
    # a_rel = g*(sin_theta - mu*cos_theta) - A_x*(cos_theta - mu*sin_theta)
    
    a_rel = g * (sin_theta - mu * cos_theta) - A_x * (cos_theta - mu * sin_theta)

    # If a_rel is not positive, the block won't slide down on its own.
    if a_rel <= 0:
        print("The block does not slide down as friction is too high or the angle is too low.")
        return

    # Calculate the distance 'L' the block slides down the slope
    L = h / sin_theta

    # Use the kinematic equation L = v_initial*t + 0.5*a_rel*t^2 to find the time 't'.
    # The block starts from rest relative to the wedge, so v_initial = 0.
    # The equation simplifies to t = sqrt(2 * L / a_rel).
    time_to_slide = math.sqrt(2 * L / a_rel)

    # Print the final calculation steps as requested
    print("The final time is calculated using the kinematic equation: L = (1/2) * a_rel * t^2")
    print("\nFirst, we find the necessary values:")
    print(f"Distance to slide, L = h / sin(theta) = {h:.2f} m / {sin_theta:.3f} = {L:.4f} m")
    print(f"Relative acceleration down the slope, a_rel = {a_rel:.4f} m/s^2")
    
    print("\nNow, we rearrange the kinematic equation to solve for t:")
    print(f"t = sqrt(2 * L / a_rel)")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"t = sqrt({(2 * L):.4f} / {a_rel:.4f})")
    print(f"t = sqrt({(2 * L / a_rel):.4f})")
    print(f"The time for the block to slide down is {time_to_slide:.4f} seconds.")
    
    # Output the final answer in the required format
    # The question asks for the exact amount of time, but the expression is
    # very complex. A high-precision float is the practical answer.
    print(f"\n<<<{time_to_slide}>>>")

solve()