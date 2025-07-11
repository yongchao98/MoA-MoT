import math

def solve_toilet_paper_problem():
    """
    Solves the falling toilet paper problem using a numerical simulation (RK4).
    """

    # 1. Define physical constants from the problem in SI units.
    G = 9.8  # Acceleration due to gravity (m/s^2)
    M_PAPER = 0.200  # Mass of the paper (kg)
    M_CYLINDER = 0.020  # Mass of the cardboard cylinder (kg)
    D_CYLINDER = 0.04  # Diameter of the cylinder (m)
    THICKNESS = 0.0005  # Thickness of one layer of paper (m)
    N_WRAPS = 100  # Number of times the paper is wrapped

    # 2. Calculate derived parameters of the system.
    R_C = D_CYLINDER / 2.0  # Radius of the cardboard cylinder (m)
    
    # The initial outer radius is the cylinder radius plus the thickness of all layers.
    R_0 = R_C + N_WRAPS * THICKNESS
    
    # The total length of the paper is calculated by conserving volume.
    # The side-view area (L * thickness) must equal the top-view area (pi * (R_0^2 - R_C^2)).
    L = math.pi * (R_0**2 - R_C**2) / THICKNESS

    # Function to calculate the roll's acceleration at a given fallen distance 'y'.
    def get_acceleration(y, v):
        # Prevent calculations with y >= L to avoid division by zero or other errors.
        if y >= L:
            return 0.0

        # Calculate instantaneous radius r(y)
        # Based on volume conservation: area_remaining = pi*(r^2 - r_c^2)
        # and area_remaining = (L-y)*thickness
        r_sq = R_C**2 + (L - y) * THICKNESS / math.pi
        r = math.sqrt(r_sq)

        # Calculate instantaneous mass M(y) of the falling roll
        fraction_remaining = (L - y) / L
        m_paper_rem = M_PAPER * fraction_remaining
        M = M_CYLINDER + m_paper_rem

        # Calculate instantaneous moment of inertia I(y)
        # Model cylinder as a thin hoop (I=mr^2) and paper as a hollow cylinder.
        I_cyl = M_CYLINDER * R_C**2
        I_paper = 0.5 * m_paper_rem * (r_sq + R_C**2)
        I = I_cyl + I_paper

        if M <= 0 or r_sq <= 0:
            return 0.0
            
        # The acceleration is g / (1 + I / (M * r^2))
        acceleration = G / (1.0 + I / (M * r_sq))
        return acceleration

    # 3. Use the 4th-order Runge-Kutta (RK4) method to solve the differential equation.
    y = 0.0  # Initial position (m)
    v = 0.0  # Initial velocity (m/s)
    t = 0.0  # Initial time (s)
    h = 0.001  # Time step for the simulation (s)

    y_prev, t_prev = 0.0, 0.0

    # Loop until the entire paper is unrolled (y >= L)
    while y < L:
        y_prev, t_prev = y, t

        # Standard RK4 integration for the system dS/dt = [v, a(y)], where S = [y, v]
        k1_v = get_acceleration(y, v)
        k1_y = v
        
        k2_v = get_acceleration(y + 0.5 * h * k1_y, v + 0.5 * h * k1_v)
        k2_y = v + 0.5 * h * k1_v
        
        k3_v = get_acceleration(y + 0.5 * h * k2_y, v + 0.5 * h * k2_v)
        k3_y = v + 0.5 * h * k2_v
        
        k4_v = get_acceleration(y + h * k3_y, v + h * k3_v)
        k4_y = v + h * k3_v
        
        # Update velocity and position
        v += (h / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        y += (h / 6.0) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
        t += h

    # 4. Interpolate to find the precise time when y == L for better accuracy.
    # final_time = t_prev + time_step * (distance_to_go / distance_in_last_step)
    final_time = t_prev + h * (L - y_prev) / (y - y_prev)

    # Output the parameters and the final answer
    print("--- Toilet Paper Roll Parameters ---")
    print(f"Cardboard cylinder diameter: {D_CYLINDER * 100} cm")
    print(f"Paper thickness: {THICKNESS * 1000} mm")
    print(f"Number of wraps: {N_WRAPS}")
    print(f"Mass of paper: {M_PAPER * 1000} grams")
    print(f"Mass of cardboard: {M_CYLINDER * 1000} grams")
    print(f"Calculated total paper length: {L:.2f} m")
    print("\n--- Result ---")
    print("The final equation for acceleration 'a' as a function of unrolled length 'y' is complex.")
    print("It is given by: a(y) = g / (1 + I(y) / (M(y) * r(y)^2))")
    print("This was solved using numerical simulation.")
    print(f"\nTime to unroll completely: {final_time:.2f} seconds")
    print(f"<<<{final_time:.2f}>>>")

solve_toilet_paper_problem()