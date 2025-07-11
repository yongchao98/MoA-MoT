import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a falling toilet paper roll to completely unroll.
    The calculation is done by numerically solving the equations of motion using the RK4 method.
    """
    # --- 1. Define physical constants from the problem ---
    G = 9.81  # Acceleration due to gravity (m/s^2)
    D_C = 4.0 / 100  # Cardboard cylinder diameter (m)
    R_C = D_C / 2.0  # Cardboard cylinder radius (m)
    T_PAPER = 0.5 / 1000  # Paper thickness (m)
    N_WRAPS = 100  # Number of wraps
    M_P = 200.0 / 1000  # Mass of paper (kg)
    M_C = 20.0 / 1000  # Mass of cardboard cylinder (kg)

    # --- 2. Calculate derived properties of the roll ---
    # Initial outer radius of the full roll
    R_O = R_C + N_WRAPS * T_PAPER

    # Total length of the paper, calculated from its volume (Area * thickness)
    # Total Area = pi * (R_o^2 - R_c^2)
    # Total Length = Total Area / thickness
    L_TOTAL = np.pi * (R_O**2 - R_C**2) / T_PAPER

    # Moment of inertia of the cardboard cylinder (approximated as a thin hoop)
    I_C = M_C * R_C**2
    
    # Pre-calculate some values for the acceleration function to optimize
    R_O_SQ = R_O**2
    R_C_SQ = R_C**2
    # This denominator is proportional to the total area of the paper
    PAPER_AREA_DENOM = R_O_SQ - R_C_SQ

    # --- 3. Define the function for instantaneous acceleration ---
    def get_acceleration(y):
        """
        Calculates the downward acceleration 'a' of the roll given the unrolled length y.
        As y increases, the roll's mass, radius, and moment of inertia decrease.
        """
        # If all paper is unrolled, the roll is just the cardboard tube.
        # This is the limiting case as y approaches L_TOTAL.
        if y >= L_TOTAL or PAPER_AREA_DENOM <= 0:
            # Acceleration of a simple hoop unrolling
            return G / 2.0

        # Calculate current radius squared, r(y)^2, based on unrolled length y
        r_sq = R_O_SQ - (y * T_PAPER) / np.pi
        r_sq = max(r_sq, R_C_SQ) # Prevent radius from going below the core

        # Calculate mass of the paper remaining on the roll
        m_p_remaining = M_P * (r_sq - R_C_SQ) / PAPER_AREA_DENOM
        m_p_remaining = max(0.0, m_p_remaining)

        # Total instantaneous mass of the roll
        M_total = M_C + m_p_remaining

        # Moment of inertia of the remaining paper (a hollow cylinder)
        I_p_remaining = 0.5 * m_p_remaining * (r_sq + R_C_SQ)
        
        # Total instantaneous moment of inertia
        I_total = I_C + I_p_remaining

        # The core equation of motion: a = g / (1 + I / (M * r^2))
        if M_total <= 0 or r_sq <= 0:
             return G # Should not happen, but as a fallback, freefall

        acceleration = G / (1 + I_total / (M_total * r_sq))
        return acceleration

    # --- 4. Run the RK4 simulation ---
    # Initial conditions
    y = 0.0  # distance fallen (m)
    v = 0.0  # velocity (m/s)
    t = 0.0  # time (s)
    dt = 0.0005  # time step for simulation (s)

    while y < L_TOTAL:
        # Standard RK4 method to solve for y(t) and v(t)
        # from dv/dt = a(y) and dy/dt = v
        k1_v = v
        k1_a = get_acceleration(y)
        
        k2_v = v + 0.5 * dt * k1_a
        k2_a = get_acceleration(y + 0.5 * dt * k1_v)
        
        k3_v = v + 0.5 * dt * k2_a
        k3_a = get_acceleration(y + 0.5 * dt * k2_v)
        
        k4_v = v + dt * k3_a
        k4_a = get_acceleration(y + dt * k3_v)

        # Update position and velocity
        y += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        v += (dt / 6.0) * (k1_a + 2*k2_a + 2*k3_a + k4_a)
        
        t += dt

    # The loop stops when y > L_TOTAL. Correct for the slight time overshoot.
    overshoot_dist = y - L_TOTAL
    time_adjustment = overshoot_dist / v if v > 0 else 0
    t -= time_adjustment
    
    # --- 5. Print the results ---
    print("Problem Parameters and Calculated Values:")
    print(f"  - Cardboard cylinder diameter: {D_C * 100:.2f} cm")
    print(f"  - Paper thickness: {T_PAPER * 1000:.2f} mm")
    print(f"  - Number of wraps: {N_WRAPS}")
    print(f"  - Mass of paper: {M_P * 1000:.0f} g, Mass of cylinder: {M_C * 1000:.0f} g")
    print(f"  - Total length of paper: {L_TOTAL:.2f} m")
    print(f"  - Initial acceleration: {get_acceleration(0):.2f} m/s^2")
    print(f"  - Final acceleration: {get_acceleration(L_TOTAL):.2f} m/s^2")
    print("\nFinal Result:")
    # Fulfilling the "output each number in the final equation" request
    print(f"Time to unroll = {t:.2f} s")

# Execute the solver function
solve_toilet_paper_fall()