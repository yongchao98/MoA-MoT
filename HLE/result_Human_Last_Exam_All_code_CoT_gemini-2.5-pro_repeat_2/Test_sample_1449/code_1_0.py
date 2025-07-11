import math

def solve():
    """
    This function calculates the time it takes for a block to slide down a moving wedge.
    """
    # Given parameters
    m = 0.1  # kg
    M = 10.0  # kg
    theta_deg = 30.0  # degrees
    h = 2.0  # m
    mu = 0.5  # coefficient of friction
    g = 10.0  # m/s^2
    v0 = 1.0 # m/s (This initial velocity does not affect the acceleration)

    # Convert angle to radians for math functions
    theta = math.radians(theta_deg)
    sin_th = math.sin(theta)
    cos_th = math.cos(theta)

    # Step 1: Solve for the acceleration of the wedge (a_M) using the non-inertial frame analysis.
    # The formula for a_M is derived from the force balance on the block and the wedge.
    # a_M = (m*g*cos(theta)*(sin(theta) + mu*cos(theta))) / (M + m*sin(theta)**2 + m*mu*sin(theta)*cos(theta))
    a_M_numerator = m * g * cos_th * (sin_th + mu * cos_th)
    a_M_denominator = M + m * sin_th**2 + m * mu * sin_th * cos_th
    a_M = a_M_numerator / a_M_denominator

    # Step 2: Solve for the relative acceleration of the block down the wedge (a_rel).
    # The formula for a_rel is also derived from the force balance in the non-inertial frame.
    # a_rel = g*(sin(theta) - mu*cos(theta)) + a_M*(mu*sin(theta) - cos(theta))
    term1 = g * (sin_th - mu * cos_th)
    term2 = a_M * (mu * sin_th - cos_th)
    a_rel = term1 + term2

    # Step 3: Calculate the length of the slope (L).
    L = h / sin_th

    # Step 4: Calculate the time (t) using kinematics.
    # L = v_initial_rel * t + 0.5 * a_rel * t^2
    # Since v_initial_rel = 0, t = sqrt(2 * L / a_rel)
    time = math.sqrt(2 * L / a_rel)

    # Output the results, including the final equation with numbers.
    print(f"The length of the slope the block slides down is L = {h} / sin({theta_deg}) = {L:.4f} m.")
    print(f"The acceleration of the wedge is a_M = {a_M:.4f} m/s^2.")
    print(f"The relative acceleration of the block down the slope is a = {a_rel:.4f} m/s^2.")
    print("\nThe time 't' is found using the kinematic equation: L = (1/2) * a * t^2")
    print("Solving for t: t = sqrt(2 * L / a)")
    print(f"Plugging in the values: t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"The time it takes for the block to slide to the bottom is t = {time:.4f} s.")
    
    # Returning the final numerical answer as a float
    return time

# Execute the function and print the final answer in the required format.
final_time = solve()
print(f"<<<{final_time:.4f}>>>")