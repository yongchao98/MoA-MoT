import numpy as np

def solve_toilet_paper_problem():
    """
    Solves the falling toilet paper problem using RK4 numerical integration.
    """
    # --- Problem Parameters ---
    g = 9.81  # Acceleration due to gravity (m/s^2)
    d_inner = 0.04  # Diameter of cardboard cylinder (m)
    paper_thickness = 0.0005  # Thickness of the paper (m)
    n_wraps = 100  # Number of wraps
    m_paper_total = 0.2  # Total mass of the paper (kg)
    m_cardboard = 0.02  # Mass of the cardboard cylinder (kg)

    # --- Derived Constants ---
    r_inner = d_inner / 2.0
    r_outer = r_inner + n_wraps * paper_thickness
    
    # Total length of the paper using a continuous model
    # Area = pi*(r_outer^2 - r_inner^2), Length = Area / thickness
    L_total = np.pi * (r_outer**2 - r_inner**2) / paper_thickness
    
    # Cross-sectional area of the paper
    A_paper_cross_section = np.pi * (r_outer**2 - r_inner**2)
    
    # Mass density per unit cross-sectional area (kg/m^2)
    sigma = m_paper_total / A_paper_cross_section
    
    # Moment of inertia of the cardboard cylinder (treated as a thin hoop)
    I_cardboard = m_cardboard * r_inner**2

    print("--- System Parameters ---")
    print(f"Total paper length (L_total): {L_total:.4f} m")
    print(f"Outer radius of full roll (r_outer): {r_outer:.4f} m")
    print(f"Inner radius of roll (r_inner): {r_inner:.4f} m")
    print(f"Moment of inertia of cardboard (I_cardboard): {I_cardboard:.2e} kg*m^2")
    print("-" * 25)

    def derivatives(state):
        """
        Calculates the derivatives [dy/dt, dv/dt] for the RK4 solver.
        state is a numpy array [y, v].
        """
        y, v = state

        # If fully unrolled, acceleration is zero (loop condition should prevent this)
        if y >= L_total:
            return np.array([v, 0.0])

        # Calculate current properties of the roll
        # Current outer radius of the paper squared
        r_sq = r_outer**2 - y * paper_thickness / np.pi
        # Clamp to r_inner^2 to avoid numerical errors at the end
        r_sq = max(r_sq, r_inner**2)
        r = np.sqrt(r_sq)

        # Mass of remaining paper
        m_paper_rem = sigma * np.pi * (r_sq - r_inner**2)
        m_paper_rem = max(0, m_paper_rem)
        
        # Total mass of the falling object
        M_total = m_cardboard + m_paper_rem

        # Moment of inertia of remaining paper (hollow cylinder)
        I_paper_rem = 0.5 * sigma * np.pi * (r_sq**2 - r_inner**4)
        I_paper_rem = max(0, I_paper_rem)
        
        # Total moment of inertia
        I_total = I_cardboard + I_paper_rem

        # Acceleration a(y, v)
        # a = (M*g - T) / M
        # T*r = I*alpha
        # a = d(r*w)/dt leads to a = r*alpha - t*v^2/(2*pi*r^2)
        # Solving for a:
        numerator = M_total * g - (I_total * paper_thickness * v**2) / (2 * np.pi * r**4)
        denominator = M_total + I_total / r**2
        
        if denominator < 1e-9: # Avoid division by zero
            a = 0
        else:
            a = numerator / denominator

        return np.array([v, a])

    # --- RK4 Solver ---
    t = 0.0
    state = np.array([0.0, 0.0])  # Initial state [y, v]
    dt = 0.001  # Time step (s)

    t_prev, y_prev = 0.0, 0.0

    while state[0] < L_total:
        t_prev, y_prev = t, state[0]

        k1 = dt * derivatives(state)
        k2 = dt * derivatives(state + 0.5 * k1)
        k3 = dt * derivatives(state + 0.5 * k2)
        k4 = dt * derivatives(state + k3)
        
        state = state + (k1 + 2*k2 + 2*k3 + k4) / 6.0
        t += dt

    # --- Final Time Calculation ---
    # Interpolate to find the exact time when y = L_total
    y_curr = state[0]
    t_curr = t
    
    if (y_curr - y_prev) > 1e-9:
        t_final = t_prev + (t_curr - t_prev) * (L_total - y_prev) / (y_curr - y_prev)
    else:
        t_final = t_curr
        
    print("\n--- Final Time Calculation ---")
    print("The simulation finished when the unrolled length passed the total length.")
    print("The final time is found by linear interpolation between the last two steps:")
    print(f"Final Time = t_prev + (t_curr - t_prev) * (L_total - y_prev) / (y_curr - y_prev)")
    print(f"Final Time = {t_prev:.4f} + ({t_curr:.4f} - {t_prev:.4f}) * ({L_total:.4f} - {y_prev:.4f}) / ({y_curr:.4f} - {y_prev:.4f})")
    print(f"Final Time = {t_final:.4f} s")
    
    print("\n--- Result ---")
    print(f"The time it takes for the toilet paper to completely unroll is {t_final:.2f} seconds.")
    
    return t_final

if __name__ == '__main__':
    final_time = solve_toilet_paper_problem()
    print(f"\n<<< {final_time:.2f} >>>")
