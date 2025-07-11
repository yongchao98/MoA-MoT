import math

def solve_time_to_slide():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    """
    # Given parameters
    m = 0.1  # kg
    M = 10.0 # kg
    g = 10.0 # m/s^2
    theta_deg = 30.0
    h = 2.0  # m
    mu = 0.5

    # Convert angle to radians for math functions
    theta_rad = math.radians(theta_deg)

    # Calculate sin and cos of theta
    s_t = math.sin(theta_rad)
    c_t = math.cos(theta_rad)

    # Step 1: Calculate the distance the block slides along the incline
    s = h / s_t

    # Step 2: Derive and calculate the acceleration of the wedge (A).
    # Using a non-inertial frame on the wedge, we can derive the wedge's acceleration A.
    # The block's motion down the slope pushes the wedge horizontally.
    # The formula for the magnitude of the wedge's acceleration is:
    # A = (m * g * cos(theta) * (sin(theta) - mu * cos(theta))) / (M + m * sin(theta) * (sin(theta) - mu * cos(theta)))
    
    # Let's calculate the terms in the formula for A
    term_A_num_factor = s_t - mu * c_t
    
    A_numerator = m * g * c_t * term_A_num_factor
    A_denominator = M + m * (s_t**2 - mu * s_t * c_t)
    
    A = A_numerator / A_denominator

    # Step 3: Calculate the acceleration of the block relative to the wedge (a_rel).
    # The formula for a_rel is derived from the forces in the non-inertial frame:
    # a_rel = g * (sin(theta) - mu * cos(theta)) + A * (cos(theta) + mu * sin(theta))
    
    term1 = g * (s_t - mu * c_t)
    term2 = A * (c_t + mu * s_t)
    
    a_rel = term1 + term2

    # Step 4: Use kinematics to find the time (t).
    # s = 0.5 * a_rel * t^2  =>  t = sqrt(2 * s / a_rel)
    t = math.sqrt(2 * s / a_rel)

    # Print the components of the final calculation as requested
    print("The final time is calculated using the equation: t = sqrt(2 * s / a_rel)")
    print(f"Sliding distance (s): {s} m")
    print(f"Relative acceleration (a_rel): {a_rel} m/s^2")
    print(f"Time to slide (t) = sqrt(2 * {s:.4f} / {a_rel:.4f})")
    print(f"Time to slide (t): {t} s")


solve_time_to_slide()
<<<3.439167386022744>>>