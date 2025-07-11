import math
import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to unroll completely
    when dropped, using a numerical simulation (RK4 method).
    """
    # 1. Define physical constants from the problem description
    M_PAPER = 0.2  # Mass of paper (kg)
    M_CYLINDER = 0.02  # Mass of cardboard cylinder (kg)
    D_INNER = 0.04  # Diameter of inner cylinder (m)
    R_INNER = D_INNER / 2.0 # Radius of inner cylinder (m)
    THICKNESS = 0.0005  # Thickness of the paper (m)
    N_WRAPS = 100  # Number of times paper is wrapped around the roll
    G = 9.8  # Gravitational acceleration (m/s^2)
    H_INITIAL = 50 # Initial height (m), not needed for this question

    # 2. Calculate derived constants
    # Total length L is the sum of circumferences of all 100 layers
    # L = Σ 2*π*r_n where r_n = R_INNER + n*THICKNESS
    L_TOTAL = 2 * math.pi * (N_WRAPS * R_INNER + THICKNESS * N_WRAPS * (N_WRAPS + 1) / 2)

    # Linear mass density (mass per unit length) of the paper
    RHO_L = M_PAPER / L_TOTAL

    # 3. Define the system of Ordinary Differential Equations (ODEs)
    # The state of the system is S = [y, v], where y is unrolled length and v is velocity.
    # We need a function that returns dS/dt = [v, a].
    def model(t, S):
        y, v = S[0], S[1]

        # Stop the simulation if the roll is fully unrolled
        if y >= L_TOTAL:
            return np.array([0.0, 0.0])
        
        y_capped = min(y, L_TOTAL)

        # Calculate current properties of the roll (radius, mass, moment of inertia)
        
        # Current outer radius squared, R(y)^2
        # Based on the remaining paper area: (L_total - y)*t = π*(R(y)^2 - R_inner^2)
        r_sq_current = R_INNER**2 + (L_TOTAL - y_capped) * THICKNESS / math.pi
        
        # Mass of remaining paper
        m_paper_rem = M_PAPER * (L_TOTAL - y_capped) / L_TOTAL
        
        # Total mass of the roll (cylinder + remaining paper)
        m_roll = M_CYLINDER + m_paper_rem

        # Moment of inertia of the roll
        # I_cylinder is treated as a thin hoop: M*R^2
        i_cyl = M_CYLINDER * R_INNER**2
        # I_paper is a thick hollow cylinder: 1/2*M*(R_outer^2 + R_inner^2)
        i_paper_rem = 0.5 * m_paper_rem * (r_sq_current + R_INNER**2)
        i_roll = i_cyl + i_paper_rem

        # Calculate acceleration a = dv/dt
        # a = (M_r*g + rho_L*v^2) / (M_r + I/R^2)
        denominator = m_roll + i_roll / r_sq_current
        numerator = m_roll * G + RHO_L * v**2
        
        a = numerator / denominator

        return np.array([v, a])

    # 4. Solve the ODE using the RK4 method
    t = 0.0
    S = np.array([0.0, 0.0])  # Initial state: y(0)=0, v(0)=0
    dt = 0.001  # Time step for simulation (s)

    # Store previous state for final interpolation
    S_prev = S
    t_prev = t

    while S[0] < L_TOTAL:
        S_prev, t_prev = S, t
        # RK4 step
        k1 = model(t, S)
        k2 = model(t + 0.5 * dt, S + 0.5 * dt * k1)
        k3 = model(t + 0.5 * dt, S + 0.5 * dt * k2)
        k4 = model(t + dt, S + dt * k3)
        S = S + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t = t + dt

    # 5. Interpolate for precise final time
    # We have the last point (t_prev, S_prev[0]) where y < L_TOTAL
    # and the first point (t, S[0]) where y > L_TOTAL
    # Use linear interpolation to find the time t_final when y = L_TOTAL
    y_prev, y_curr = S_prev[0], S[0]
    t_curr = t
    
    final_time = t_prev + (t_curr - t_prev) * (L_TOTAL - y_prev) / (y_curr - y_prev)
    
    print(f"Total length of paper to unroll: {L_TOTAL:.2f} m")
    print(f"Time to unroll completely: {final_time:.2f} s")
    
    return round(final_time, 2)

final_answer = solve_toilet_paper_fall()
print(f"\nFinal Answer: The time it takes for the toilet paper to reach the end of its roll is {final_answer} seconds.")
# The final answer is wrapped in <<<>>> as requested.
print(f"<<<{final_answer}>>>")
