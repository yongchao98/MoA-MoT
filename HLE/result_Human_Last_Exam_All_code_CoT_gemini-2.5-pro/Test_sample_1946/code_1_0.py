import numpy as np

def solve_toilet_paper_problem():
    """
    Solves the falling toilet paper problem using a Runge-Kutta (RK4) numerical method.
    """
    # --- Problem Parameters ---
    g = 9.81  # Gravitational acceleration (m/s^2)
    d_c = 4.0  # Diameter of cardboard cylinder (cm)
    m_c_grams = 20.0  # Mass of cardboard cylinder (grams)
    t_paper_mm = 0.5  # Thickness of paper (mm)
    N_wraps = 100  # Number of wraps
    m_p_total_grams = 200.0  # Total mass of paper (grams)

    # --- Convert to SI Units ---
    r_c = (d_c / 2) / 100  # m (cardboard inner radius)
    m_c = m_c_grams / 1000  # kg (cardboard mass)
    t_paper = t_paper_mm / 1000  # m (paper thickness)
    m_p_total = m_p_total_grams / 1000  # kg (total paper mass)

    # --- Derived Constants ---
    # Initial outer radius of the full roll
    r_initial = r_c + N_wraps * t_paper

    # Total length of the paper, derived by modeling the paper as a spiral
    # Formula: L = (pi/t) * (r_outer^2 - r_inner^2)
    L_total = (np.pi / t_paper) * (r_initial**2 - r_c**2)

    # Mass per unit length of the paper
    lambda_p = m_p_total / L_total

    # Moment of inertia of the cardboard cylinder (modeled as a thin-walled cylinder: I = mr^2)
    I_c = m_c * r_c**2

    # --- State-dependent acceleration function ---
    def get_acceleration(y):
        """Calculates instantaneous acceleration based on the length of paper unrolled (y)."""
        # Ensure we don't calculate past the end of the roll
        if y >= L_total:
            y = L_total
        
        # Square of the current outer radius as a function of unrolled length y
        r_sq = r_initial**2 - (t_paper / np.pi) * y
        
        # If rounding errors push r_sq below r_c^2, clamp it.
        if r_sq < r_c**2:
            r_sq = r_c**2

        # Current mass of the remaining paper
        m_p_current = lambda_p * (L_total - y)

        # Current total mass of the roll
        m_current = m_c + m_p_current

        # Moment of inertia of the remaining paper (modeled as a hollow cylinder)
        # I_p = 0.5 * m * (r_outer^2 + r_inner^2)
        I_p_current = 0.5 * m_p_current * (r_sq + r_c**2)

        # Total current moment of inertia
        I_current = I_c + I_p_current

        # Calculate acceleration using the formula: a = g / (1 + I / (m * r^2))
        if m_current * r_sq == 0: # Avoid division by zero at the very end
             return g / (1 + I_c / (m_c * r_c**2))
        
        acceleration = g / (1 + I_current / (m_current * r_sq))
        
        return acceleration

    # --- RK4 Simulation ---
    y = 0.0  # Initial unrolled length
    v = 0.0  # Initial velocity
    time = 0.0
    dt = 0.001  # Time step in seconds

    # Store previous state for final interpolation
    y_prev, time_prev = y, time

    while y < L_total:
        y_prev, time_prev = y, time
        
        # RK4 steps for the ODE system: dS/dt = [v, a(y)] where S = [y, v]
        k1_v = get_acceleration(y)
        k1_y = v

        k2_v = get_acceleration(y + 0.5 * dt * k1_y)
        k2_y = v + 0.5 * dt * k1_v

        k3_v = get_acceleration(y + 0.5 * dt * k2_y)
        k3_y = v + 0.5 * dt * k2_v

        k4_v = get_acceleration(y + dt * k3_y)
        k4_y = v + dt * k3_v
        
        # Update velocity and position
        v += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        y += (dt / 6.0) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
        time += dt

    # --- Final Interpolation and Output ---
    # Interpolate to find the exact time when y = L_total
    overshoot = y - L_total
    last_step_dy = y - y_prev
    fraction = overshoot / last_step_dy
    final_time = time - fraction * dt
    
    # Calculate the values for the final acceleration equation
    m_final = m_c
    r_sq_final = r_c**2
    I_final = I_c
    final_acceleration = g / (1 + I_final / (m_final * r_sq_final))

    print(f"Total length of the toilet paper is {L_total:.2f} m.")
    print("The final acceleration is determined by the cardboard cylinder alone.")
    print("Equation: a = g / (1 + I_final / (m_final * r_final^2))")
    print(f"Values: a = {g} / (1 + {I_final:.8f} / ({m_final} * {r_sq_final:.4f})) = {final_acceleration:.2f} m/s^2")
    print("-" * 50)
    print(f"The time it takes for the toilet paper to reach the end of its roll is {final_time:.2f} s.")

solve_toilet_paper_problem()
<<<3.21>>>