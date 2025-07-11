import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time for a falling toilet paper roll to unravel using the RK4 method.
    """

    # --- 1. Physical Constants and Parameters ---
    G = 9.8  # Gravitational acceleration (m/s^2)
    D_INNER_CM = 4.0  # Inner cylinder diameter (cm)
    T_PAPER_MM = 0.5  # Paper thickness (mm)
    N_WRAPS = 100  # Number of wraps
    M_P_G = 200.0  # Total mass of paper (g)
    M_C_G = 20.0  # Mass of cardboard cylinder (g)

    # Convert to SI units (meters, kilograms)
    R_INNER = (D_INNER_CM / 2) / 100
    T_PAPER_M = T_PAPER_MM / 1000
    M_P_TOTAL_KG = M_P_G / 1000
    M_C_KG = M_C_G / 1000
    
    # --- 2. Derived Properties of the Roll ---
    # Initial outer radius of the full roll
    R_INITIAL = R_INNER + N_WRAPS * T_PAPER_M
    
    # Total length of the paper, using a continuous model for self-consistency
    # L_total = Area / thickness
    A_paper_cross_section = np.pi * (R_INITIAL**2 - R_INNER**2)
    L_TOTAL = A_paper_cross_section / T_PAPER_M
    
    # --- 3. Derivative function for the RK4 solver ---
    # Pre-calculate constant terms for efficiency
    R_INITIAL_SQ = R_INITIAL**2
    R_INNER_SQ = R_INNER**2
    DENOM_M_P = R_INITIAL_SQ - R_INNER_SQ

    def get_derivatives(state):
        """
        Calculates the derivatives [dy/dt, dv/dt] for a given state [y, v].
        
        Args:
            state: A numpy array [y, v], where y is the unrolled distance
                   and v is the downward velocity.
                   
        Returns:
            A numpy array containing the derivatives [v, a].
        """
        y, v = state
        
        # Calculate current properties of the roll based on unrolled length y
        # R(y)^2 = R_initial^2 - y*t/pi (from L = Area/t => y = pi*(Ri^2-R(y)^2)/t)
        current_R_sq = R_INITIAL_SQ - y * T_PAPER_M / np.pi
        
        # Guard against numerical errors when y is near L_TOTAL
        if current_R_sq < R_INNER_SQ:
            current_R_sq = R_INNER_SQ

        # Mass of paper remaining (proportional to remaining area)
        m_p_remaining = M_P_TOTAL_KG * (current_R_sq - R_INNER_SQ) / DENOM_M_P
        
        # Total instantaneous mass of the falling roll
        M_total = M_C_KG + m_p_remaining
        
        # Moment of inertia
        # Cardboard: modeled as a thin ring, I_c = M_c * r^2
        I_c = M_C_KG * R_INNER_SQ
        # Paper: modeled as a hollow cylinder, I_p = 1/2 * m_p * (R_outer^2 + R_inner^2)
        I_p = 0.5 * m_p_remaining * (current_R_sq + R_INNER_SQ)
        I_total = I_c + I_p
        
        # Acceleration a = g / (1 + I / (M * R^2))
        # If the roll is a simple object, the denominator is a constant (e.g. 1.5 for a solid cylinder).
        # Here, it changes as the roll unrolls.
        a = G / (1 + I_total / (M_total * current_R_sq))
        
        return np.array([v, a])

    # --- 4. RK4 Solver Implementation ---
    # Initial conditions
    t = 0.0
    # State vector S = [y, v] = [position, velocity]
    state = np.array([0.0, 0.0])
    dt = 0.001  # Time step in seconds

    # Store previous state for final interpolation
    t_prev, state_prev = t, np.copy(state)

    while state[0] < L_TOTAL:
        # Store state before the step
        t_prev, state_prev = t, np.copy(state)
        
        # RK4 calculation steps
        k1 = get_derivatives(state)
        k2 = get_derivatives(state + 0.5 * dt * k1)
        k3 = get_derivatives(state + 0.5 * dt * k2)
        k4 = get_derivatives(state + dt * k3)
        
        # Update state and time
        state = state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t = t + dt

    # --- 5. Final Result ---
    # Interpolate to find the exact time when y = L_TOTAL
    y_prev = state_prev[0]
    y_curr = state[0]
    # t_final = t_prev + dt * (L_TOTAL - y_prev) / (y_curr - y_prev)
    t_final = t_prev + (L_TOTAL - y_prev) * (t - t_prev) / (y_curr - y_prev)

    # --- 6. Output the numbers and the result ---
    print("--- Input Parameters ---")
    print(f"Cardboard Mass (M_c): {M_C_KG:.3f} kg")
    print(f"Total Paper Mass (M_p): {M_P_TOTAL_KG:.3f} kg")
    print(f"Inner Radius (r_in): {R_INNER:.3f} m")
    print(f"Paper Thickness (t): {T_PAPER_M} m")
    print("\n--- Calculated Properties ---")
    print(f"Initial Outer Radius (R_initial): {R_INITIAL:.3f} m")
    print(f"Total Length of Paper (L_total): {L_TOTAL:.2f} m")
    
    a_start = get_derivatives(np.array([0.0, 0.0]))[1]
    a_end = get_derivatives(np.array([L_TOTAL, 0.0]))[1]
    print(f"Initial Acceleration: {a_start:.2f} m/s^2")
    print(f"Final Acceleration: {a_end:.2f} m/s^2")

    print("\n--- Final Answer ---")
    print(f"The time it takes for the toilet paper to completely unroll is {t_final:.2f} seconds.")
    
# Execute the solver function
solve_toilet_paper_fall()
<<<3.14>>>