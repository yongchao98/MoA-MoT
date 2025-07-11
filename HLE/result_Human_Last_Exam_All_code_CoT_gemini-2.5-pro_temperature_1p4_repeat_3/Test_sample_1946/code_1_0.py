import numpy as np
from scipy.integrate import solve_ivp

def solve_toilet_paper_problem():
    """
    Calculates the time it takes for a falling toilet paper roll to completely unroll.

    The solution is found by numerically solving the equations of motion derived
    from the principle of conservation of energy using a Runge-Kutta method.
    """
    # 1. Define Physical Parameters
    # All units are in SI (meters, kilograms, seconds).
    G = 9.80665              # Acceleration due to gravity (m/s^2)
    R_INNER = 0.04 / 2       # Inner radius of cardboard cylinder (m)
    THICKNESS = 0.5 / 1000   # Thickness of the paper (m)
    N_WRAPS = 100            # Number of wraps
    M_PAPER = 200 / 1000     # Mass of the paper (kg)
    M_CYLINDER = 20 / 1000   # Mass of the cardboard cylinder (kg)

    # Derived constants
    R_OUTER = R_INNER + N_WRAPS * THICKNESS
    # Total length of the paper, calculated from the cross-sectional area
    L_TOTAL = np.pi * (R_OUTER**2 - R_INNER**2) / THICKNESS
    # Linear mass density of the paper
    LAMBDA_L = M_PAPER / L_TOTAL
    # Moment of inertia of the cylinder (modeled as a thin hoop, I = MR^2)
    I_CYLINDER = M_CYLINDER * R_INNER**2

    # Using a dictionary for memoization to speed up the ODE solver
    memo = {}

    def get_acceleration(y):
        """
        Calculates the acceleration of the roll as a function of unrolled length y.
        This function is derived from the conservation of energy principle.
        a = 0.5 * d(v^2)/dy
        """
        if y in memo:
            return memo[y]
        
        # Handle the initial state to avoid division by zero
        if y < 1e-9:
            i_0 = I_CYLINDER + 0.5 * M_PAPER * (R_INNER**2 + R_OUTER**2)
            m_0 = M_CYLINDER + M_PAPER
            return (m_0 * G) / (m_0 + i_0 / R_OUTER**2)

        # Calculate instantaneous properties of the roll at unrolled length y
        m_paper_rem_y = LAMBDA_L * (L_TOTAL - y)
        m_y = M_CYLINDER + m_paper_rem_y
        r_sq_y = R_OUTER**2 - y * THICKNESS / np.pi
        
        # Ensure radius squared is positive
        if r_sq_y < 0: r_sq_y = 1e-12 

        # Instantaneous moment of inertia of the paper (as a hollow cylinder)
        i_paper_y = 0.5 * m_paper_rem_y * (R_INNER**2 + r_sq_y)
        i_y = I_CYLINDER + i_paper_y

        # v^2 = N(y) / D(y)
        n_y = 2 * G * (m_y * y + 0.5 * LAMBDA_L * y**2)
        d_y = m_y + i_y / r_sq_y + LAMBDA_L * y

        if d_y == 0: return 0

        # Derivatives w.r.t. y needed for a = 0.5 * d(v^2)/dy
        dn_dy = 2 * G * m_y

        dr_sq_dy = -THICKNESS / np.pi
        di_dy = 0.5 * (-LAMBDA_L * (R_INNER**2 + r_sq_y) + m_paper_rem_y * dr_sq_dy)
        dd_dy = (di_dy * r_sq_y - i_y * dr_sq_dy) / (r_sq_y**2)
        
        # Final acceleration formula: a = 0.5 * (N'*D - N*D') / D^2
        acceleration = 0.5 * (dn_dy * d_y - n_y * dd_dy) / (d_y**2)
        
        memo[y] = acceleration
        return acceleration

    def odes(t, state):
        """Defines the system of ordinary differential equations for the solver."""
        y, v = state
        # Velocity must be non-negative.
        v_safe = max(v, 0)
        dv_dt = get_acceleration(y)
        return [v_safe, dv_dt]

    # Event function to stop the integration when the roll is fully unrolled
    def end_of_roll(t, state):
        return state[0] - L_TOTAL
    end_of_roll.terminal = True
    end_of_roll.direction = 1

    # Initial conditions [y(0), v(0)] and time span for the simulation
    initial_state = [0.0, 0.0]
    t_span = [0, 20] # A reasonable time span in seconds

    # 3. Solve the ODE system
    solution = solve_ivp(
        odes,
        t_span,
        initial_state,
        events=end_of_roll,
        dense_output=True,
        method='RK45'
    )

    # 4. Output the results
    print("--- Toilet Paper Unrolling Simulation ---")
    print(f"Inner Cylinder Diameter: {R_INNER*2*100:.2f} cm")
    print(f"Total Paper Mass:        {M_PAPER*1000:.2f} g")
    print(f"Cardboard Cylinder Mass: {M_CYLINDER*1000:.2f} g")
    print(f"Total Paper Length:      {L_TOTAL:.2f} m")

    if solution.status == 1: # Status 1 means an event terminated the integration
        final_time = solution.t_events[0][0]
        # This is the "final equation" requested: a statement of the result.
        print("\nFinal Result:")
        print(f"The time for the toilet paper to completely unroll is {final_time:.2f} seconds.")
    else:
        print("\nThe simulation did not find a solution within the given time span.")
        
    return final_time if solution.status == 1 else None

# Execute the solution
final_time = solve_toilet_paper_problem()
print(f"<<<{final_time:.2f}>>>")
