import math

def solve_time_to_slide():
    """
    Calculates the time for a block to slide down a moving wedge with friction.
    """
    # 1. Define given physical constants
    m = 0.1      # mass of the block in kg (100 g)
    M = 10.0     # mass of the wedge in kg
    theta_deg = 30.0 # angle of the wedge in degrees
    h = 2.0      # height of the wedge in m
    mu = 0.5     # coefficient of kinetic friction
    g = 10.0     # acceleration due to gravity in m/s^2

    # 2. Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    s_t = math.sin(theta_rad)
    c_t = math.cos(theta_rad)

    # 3. Solve for the wedge's acceleration (A) and the Normal force (N)
    # The derivation involves solving a system of two equations from force analysis:
    # (1) M*A = N*(sin(theta) + mu*cos(theta))
    # (2) N = m*g*cos(theta) + m*A*sin(theta)
    # Solving for A yields:
    # A = (m*g*cos(theta)*(sin(theta)+mu*cos(theta))) / (M - m*sin(theta)*(sin(theta)+mu*cos(theta)))
    
    A_numerator = m * g * c_t * (s_t + mu * c_t)
    A_denominator = M - m * s_t * (s_t + mu * c_t)
    A = A_numerator / A_denominator

    # We can now find the Normal force N using equation (2)
    N = m * g * c_t + m * A * s_t

    # 4. Solve for the relative acceleration (a_rel) of the block down the incline
    # From force analysis along the incline in the wedge's frame:
    # m*a_rel = m*g*sin(theta) - Friction - (pseudo-force component)
    # m*a_rel = m*g*sin(theta) - mu*N - m*A*cos(theta)
    # a_rel = g*sin(theta) - (mu*N)/m - A*cos(theta)
    
    a_rel = g * s_t - (mu * N) / m - A * c_t

    # 5. Calculate the length of the incline (L)
    L = h / s_t

    # 6. Calculate the time (t) using kinematics
    # L = 0.5 * a_rel * t^2  => t = sqrt(2 * L / a_rel)
    t = math.sqrt((2 * L) / a_rel)

    # Print the breakdown of the calculation as requested
    print(f"The equation for the time is: t = sqrt((2 * L) / a_rel)")
    print("\n--- Calculated Intermediate Values ---")
    print(f"Length of incline, L = {h:.1f} m / sin({theta_deg} deg) = {L:.4f} m")
    print(f"Acceleration of wedge, A = {A:.4f} m/s^2")
    print(f"Normal force, N = {N:.4f} N")
    print(f"Relative acceleration of block, a_rel = {a_rel:.4f} m/s^2")
    print("\n--- Final Calculation ---")
    print(f"t = sqrt((2 * {L:.4f} m) / {a_rel:.4f} m/s^2)")
    print(f"t = sqrt({(2*L):.4f} / {a_rel:.4f}) s")
    print(f"t = sqrt({(2*L/a_rel):.4f}) s")
    print(f"t = {t:.4f} s")

    # Final answer in the required format
    print(f"\n<<<{t}>>>")

solve_time_to_slide()