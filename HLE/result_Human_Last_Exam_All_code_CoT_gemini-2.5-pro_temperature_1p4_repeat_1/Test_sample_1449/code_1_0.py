import math

def solve_slide_time():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    """
    # Given parameters
    m = 0.1   # Mass of the block in kg
    M = 10.0  # Mass of the wedge in kg
    g = 10.0  # Acceleration due to gravity in m/s^2
    h = 2.0   # Initial height of the block in m
    mu = 0.5  # Coefficient of friction
    theta_deg = 30.0 # Angle of the wedge in degrees

    # Convert angle to radians for use in trigonometric functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # Step 1: Calculate the horizontal acceleration of the wedge (A_x)
    # The forces from the block (normal and friction) push the wedge horizontally.
    # Numerator of the A_x expression
    num_A = m * g * cos_theta * (sin_theta + mu * cos_theta)
    # Denominator of the A_x expression
    den_A = M + m * (sin_theta**2 + mu * sin_theta * cos_theta)
    A_x = num_A / den_A

    # Step 2: Calculate the acceleration of the block relative to the wedge (a_rel)
    # This is the acceleration down the slope in the wedge's reference frame.
    # It includes the component of gravity, friction, and the inertial force.
    term1 = g * (sin_theta - mu * cos_theta)
    term2 = A_x * (cos_theta - mu * sin_theta)
    a_rel = term1 - term2

    # Step 3: Calculate the distance the block needs to slide (L)
    # This is the hypotenuse of the wedge's triangle.
    L = h / sin_theta

    # Step 4: Calculate the time (t) using kinematics
    # L = v0*t + 0.5*a_rel*t^2. Since v0 (relative) is 0, t = sqrt(2*L/a_rel).
    time = math.sqrt(2 * L / a_rel)

    # Print the final equation with the calculated numbers
    print(f"To find the time, we use the equation: t = sqrt(2 * L / a_rel)")
    print(f"The distance to slide is L = {L:.4f} m.")
    print(f"The relative acceleration of the block is a_rel = {a_rel:.4f} m/s^2.")
    print(f"Plugging the numbers into the equation:")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f}) = {time:.4f} s")

    # Output the final answer in the specified format
    print(f"\n<<<{time:.4f}>>>")

solve_slide_time()