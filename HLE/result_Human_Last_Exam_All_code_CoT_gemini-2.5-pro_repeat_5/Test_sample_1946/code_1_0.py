import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a falling toilet paper roll to completely unroll.
    
    The problem is modeled using a system of ordinary differential equations, which is
    solved numerically using the fourth-order Runge-Kutta (RK4) method.
    """
    
    # 1. Define constants and initial parameters
    G = 9.81  # Acceleration due to gravity (m/s^2)
    
    # Cardboard cylinder properties
    M_C = 0.020  # Mass of cardboard cylinder (kg)
    D_C = 0.04   # Diameter of cardboard cylinder (m)
    R_C = D_C / 2.0 # Radius of cardboard cylinder (m)
    
    # Toilet paper properties
    M_P_TOTAL = 0.200 # Total mass of the paper (kg)
    THICKNESS = 0.0005 # Thickness of the paper (m, 0.5 mm)
    WRAPS = 100 # Number of times the paper is wrapped around the roll
    
    # 2. Calculate derived physical constants
    # Initial full radius of the roll
    R_F = R_C + WRAPS * THICKNESS
    # Total length of the paper, calculated by equating the cross-sectional area
    # of the paper roll to the area of the unrolled paper (L * thickness).
    L_TOTAL = np.pi * (R_F**2 - R_C**2) / THICKNESS

    # 3. Define the function for the system of ODEs (for RK4)
    # The state vector is [y, v], where y is the distance fallen and v is the velocity.
    # The function returns the derivatives [dy/dt, dv/dt] = [v, a].
    def get_state_derivatives(state):
        y, v = state
        
        # Stop calculation if the roll is fully unrolled
        if y >= L_TOTAL:
            return np.array([v, 0.0])

        # Fraction of paper remaining
        frac_rem = max(0, 1.0 - y / L_TOTAL)
        
        # Calculate instantaneous mass
        m_p_rem = M_P_TOTAL * frac_rem
        m_total = M_C + m_p_rem

        # Calculate instantaneous radius squared
        # R(y)^2 = R_c^2 + (R_f^2 - R_c^2) * fraction_remaining
        r_sq = R_C**2 + (R_F**2 - R_C**2) * frac_rem
        # Prevent division by zero at the very end of the roll
        if r_sq < 1e-9:
            r_sq = 1e-9

        # Calculate instantaneous moment of inertia for the roll
        # I_cylinder = M_c * R_c^2 (thin shell approximation)
        i_c = M_C * R_C**2
        # I_paper = 1/2 * M_paper * (R_outer^2 + R_inner^2)
        i_p = 0.5 * m_p_rem * (r_sq + R_C**2)
        i_total = i_c + i_p

        # Calculate instantaneous acceleration using the formula: a = m*g / (m + I/R^2)
        acceleration = (m_total * G) / (m_total + i_total / r_sq)
        
        return np.array([v, acceleration])

    # 4. RK4 numerical integration
    t = 0.0
    dt = 0.001  # Time step for the simulation (s)
    state = np.array([0.0, 0.0])  # Initial state: y=0, v=0
    
    # Store the previous state to use for final interpolation
    t_prev, state_prev = t, state

    while state[0] < L_TOTAL:
        t_prev, state_prev = t, state
        
        k1 = get_state_derivatives(state)
        k2 = get_state_derivatives(state + 0.5 * dt * k1)
        k3 = get_state_derivatives(state + 0.5 * dt * k2)
        k4 = get_state_derivatives(state + dt * k3)
        
        state = state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t = t + dt

    # 5. Interpolate to find the exact time when y = L_TOTAL
    y_prev = state_prev[0]
    y_curr = state[0]
    t_curr = t
    
    # Avoid division by zero if the step size is too large
    if (y_curr - y_prev) > 1e-9:
        # t_final = t_prev + (t_curr - t_prev) * (L_TOTAL - y_prev) / (y_curr - y_prev)
        time_adjustment_factor = (L_TOTAL - y_prev) / (y_curr - y_prev)
        t_final = t_prev + (t_curr - t_prev) * time_adjustment_factor
    else:
        t_final = t_curr

    # 6. Output the results
    print("The final time is found by interpolating the last simulation step:")
    print(f"t_final = {t_prev:.4f} + ({t_curr:.4f} - {t_prev:.4f}) * ({L_TOTAL:.4f} - {y_prev:.4f}) / ({y_curr:.4f} - {y_prev:.4f})")
    print(f"\nTime to reach the end of the roll: {t_final:.2f} s")

solve_toilet_paper_fall()
<<<3.08>>>