import numpy as np
from scipy.integrate import solve_ivp

def solve_toilet_paper_problem():
    """
    Solves the falling toilet paper problem using numerical integration.
    """
    # --- Input Parameters ---
    # These are the numbers that define the problem.
    G = 9.81  # Acceleration due to gravity (m/s^2)
    R_C = 0.04 / 2  # Radius of the cardboard cylinder (m)
    M_C = 0.020  # Mass of the cardboard cylinder (kg)
    M_P = 0.200  # Total mass of the paper (kg)
    THICKNESS = 0.0005  # Thickness of the paper (m)
    N_WRAPS = 100  # Number of wraps
    H_INITIAL = 50 # Initial height (m)
    
    # --- Derived and Helper Constants ---
    R_INITIAL = R_C + N_WRAPS * THICKNESS
    R_C_SQR = R_C**2
    R_INITIAL_SQR = R_INITIAL**2
    R_DIFF_SQR = R_INITIAL_SQR - R_C_SQR # Denominator for mass ratio, pre-calculated
    MC_RC_SQR = M_C * R_C_SQR # I_cardboard, pre-calculated

    def derivatives(t, S):
        """
        Calculates the derivatives dS/dt for the ODE solver.
        S is the state vector [v, omega, r].
        v: vertical velocity (m/s)
        omega: angular velocity (rad/s)
        r: current outer radius (m)
        """
        _v, omega, r = S
        
        # Stop if the paper is fully unrolled to prevent invalid calculations (r < R_C)
        if r <= R_C:
            return [0, 0, 0]

        r_sqr = r**2

        # Calculate current mass m(r)
        m_paper_rem = M_P * (r_sqr - R_C_SQR) / R_DIFF_SQR
        m_current = M_C + m_paper_rem

        # Calculate current moment of inertia I(r)
        # I = I_cardboard + I_paper = (M_c * R_c^2) + (1/2 * m_paper * (r^2 + R_c^2))
        I_paper_rem = 0.5 * m_paper_rem * (r_sqr + R_C_SQR)
        I_current = MC_RC_SQR + I_paper_rem
        
        # Acceleration a = m*g / (m + I/r^2)
        acceleration = (m_current * G) / (m_current + I_current / r_sqr)
        
        # The derivatives of our state vector [v, omega, r]
        dv_dt = acceleration
        domega_dt = acceleration / r
        dr_dt = - (THICKNESS / (2 * np.pi)) * omega
        
        return [dv_dt, domega_dt, dr_dt]

    # --- Numerical Solution ---
    # Initial conditions: S = [v, omega, r] at t=0
    S0 = [0.0, 0.0, R_INITIAL]

    # Event function: triggers when r(t) = R_C
    def unrolled_event(t, S):
        return S[2] - R_C # Returns zero when r = R_C
    
    unrolled_event.terminal = True  # Stop the integration when the event occurs
    unrolled_event.direction = -1   # Trigger only when r is decreasing

    # Set a time span for the simulation, it will stop early if event is found
    t_span = (0, 20)

    # Call the ODE solver
    solution = solve_ivp(
        fun=derivatives,
        t_span=t_span,
        y0=S0,
        events=unrolled_event,
        method='RK45'
    )
    
    # --- Output Final Results ---
    print("This script solves for the time it takes a roll of toilet paper to unravel while falling.")
    print("The final calculation is based on the following equation components and parameters:")
    print(f"  Gravity (g): {G} m/s^2")
    print(f"  Cardboard Cylinder Diameter: {2 * R_C * 100:.1f} cm")
    print(f"  Mass of Cardboard Cylinder (m_c): {M_C * 1000:.0f} g")
    print(f"  Total Mass of Paper (m_p): {M_P * 1000:.0f} g")
    print(f"  Paper Thickness: {THICKNESS * 1000:.1f} mm")
    print(f"  Number of Wraps: {N_WRAPS}")

    print("\nThe problem is modeled by a system of differential equations and solved numerically.")
    print("System State S = [velocity, angular velocity, radius]")
    print("d_S/dt = [ a , a/r , -(thickness/2pi)*omega ]")
    print("where acceleration 'a' is a function of the changing mass and moment of inertia.")
    
    if solution.status == 1: # Status 1 means a terminal event was triggered
        time_to_unroll = solution.t_events[0][0]
        print(f"\nThe time for the toilet paper to completely unroll is {time_to_unroll:.2f} seconds.")
    else:
        print("\nThe solver did not find the solution within the given time span.")

solve_toilet_paper_problem()