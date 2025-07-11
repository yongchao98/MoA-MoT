import math

def solve_toilet_paper_problem():
    """
    Solves the falling toilet paper problem using a numerical simulation (RK4).
    """

    # 1. Define given constants in SI units
    D_INNER_CM = 4.0
    PAPER_THICKNESS_MM = 0.5
    N_WRAPS = 100
    M_PAPER_G = 200.0
    M_CORE_G = 20.0
    G = 9.8  # Acceleration due to gravity in m/s^2

    # Convert to SI units (meters, kilograms)
    R_INNER = (D_INNER_CM / 100) / 2.0  # m
    PAPER_THICKNESS = PAPER_THICKNESS_MM / 1000  # m
    M_PAPER = M_PAPER_G / 1000  # kg
    M_CORE = M_CORE_G / 1000  # kg

    # 2. Derive initial parameters
    R_OUTER = R_INNER + N_WRAPS * PAPER_THICKNESS
    R_OUTER_SQ = R_OUTER**2
    R_INNER_SQ = R_INNER**2

    # Calculate total length of the paper based on volume/area consistency.
    # Area = pi * (R_outer^2 - R_inner^2) = Length * Thickness
    PAPER_TOTAL_LENGTH = math.pi * (R_OUTER_SQ - R_INNER_SQ) / PAPER_THICKNESS

    def get_acceleration(y):
        """
        Calculates the acceleration 'a' for a given unrolled length 'y'.
        This is the implementation of a = m*g / (m + I/r^2).
        """
        # Stop condition: if y exceeds total length, clamp it.
        if y >= PAPER_TOTAL_LENGTH:
            y = PAPER_TOTAL_LENGTH

        # Calculate instantaneous properties of the roll
        fraction_remaining = 1.0 - y / PAPER_TOTAL_LENGTH

        # Current mass (core + remaining paper)
        mass_paper_remaining = M_PAPER * fraction_remaining
        mass = M_CORE + mass_paper_remaining

        # Current radius squared
        # r^2 = r_inner^2 + (r_outer^2 - r_inner^2) * fraction_remaining
        r_sq = R_INNER_SQ + (R_OUTER_SQ - R_INNER_SQ) * fraction_remaining
        if r_sq < 1e-12: r_sq = 1e-12 # Prevent division by zero

        # Current moment of inertia (core + remaining paper)
        # Core is a thin shell: I_core = m * r^2
        # Paper is a hollow cylinder: I_paper = 0.5 * m * (r_outer^2 + r_inner^2)
        inertia_core = M_CORE * R_INNER_SQ
        inertia_paper = 0.5 * mass_paper_remaining * (r_sq + R_INNER_SQ)
        inertia = inertia_core + inertia_paper

        # Calculate acceleration a = m*g / (m + I/r^2)
        acceleration = (mass * G) / (mass + inertia / r_sq)
        return mass, r_sq**0.5, inertia, acceleration

    # 3. Print the breakdown of the governing equation at start and end points
    print("This problem is solved by simulating the fall numerically.")
    print("The governing equation for acceleration (a) is: a = m*g / (m + I/r^2)\n")

    # Values at the start (y=0)
    m0, r0, I0, a0 = get_acceleration(0)
    print("Equation parameters at the start (y=0):")
    print(f"  m (mass)      = {m0:.3f} kg")
    print(f"  r (radius)    = {r0:.3f} m")
    print(f"  I (inertia)   = {I0:.6f} kg*m^2")
    print(f"  Acceleration a(0) = ({m0:.3f} * {G}) / ({m0:.3f} + {I0:.6f} / {r0:.3f}^2) = {a0:.2f} m/s^2\n")

    # Values at the end (y=L)
    mL, rL, IL, aL = get_acceleration(PAPER_TOTAL_LENGTH)
    print("Equation parameters at the end (y=L):")
    print(f"  m (mass)      = {mL:.3f} kg")
    print(f"  r (radius)    = {rL:.3f} m")
    print(f"  I (inertia)   = {IL:.8f} kg*m^2")
    print(f"  Acceleration a(L) = ({mL:.3f} * {G}) / ({mL:.3f} + {IL:.8f} / {rL:.3f}^2) = {aL:.2f} m/s^2\n")

    # 4. Run the RK4 simulation
    y = 0.0  # Position (distance fallen/unrolled)
    v = 0.0  # Velocity
    t = 0.0  # Time
    dt = 0.001  # Time step for simulation

    y_prev, t_prev = y, t

    while y < PAPER_TOTAL_LENGTH:
        y_prev, t_prev = y, t

        # Standard RK4 method for y'' = f(y)
        k1_v = get_acceleration(y)[3]
        k1_y = v
        k2_v = get_acceleration(y + 0.5 * dt * k1_y)[3]
        k2_y = v + 0.5 * dt * k1_v
        k3_v = get_acceleration(y + 0.5 * dt * k2_y)[3]
        k3_y = v + 0.5 * dt * k2_v
        k4_v = get_acceleration(y + dt * k3_y)[3]
        k4_y = v + dt * k3_v

        v += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        y += (dt / 6.0) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
        t += dt

    # Interpolate for the exact time when y = PAPER_TOTAL_LENGTH
    # t_final = t_prev + dt * (L - y_prev) / (y_curr - y_prev)
    final_time = t_prev + dt * (PAPER_TOTAL_LENGTH - y_prev) / (y - y_prev)
    
    # 5. Output the final answer
    print(f"The total length of the paper is {PAPER_TOTAL_LENGTH:.2f} m.")
    print("Running simulation...")
    print("\n---------------------------------------------------------")
    print(f"The time it takes to completely unroll is: {final_time:.2f} seconds")
    print("---------------------------------------------------------")
    
    return final_time

# Execute the solver
final_time_result = solve_toilet_paper_problem()