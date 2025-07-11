import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to unroll while falling.
    
    The problem is solved by numerically integrating the equations of motion
    using the 4th-order Runge-Kutta (RK4) method, as the system's mass,
    radius, and moment of inertia change continuously.
    """
    # --- Given Parameters ---
    g = 9.81  # Acceleration due to gravity (m/s^2)
    d_inner_cm = 4.0  # Diameter of inner cylinder (cm)
    paper_thickness_mm = 0.5  # Thickness of paper (mm)
    num_wraps = 100  # Number of wraps
    mass_paper_g = 200.0  # Total mass of paper (g)
    mass_cylinder_g = 20.0  # Mass of cardboard cylinder (g)

    # --- Convert to SI Units ---
    r_i = (d_inner_cm / 100) / 2      # Inner radius (m)
    tau = paper_thickness_mm / 1000   # Paper thickness (m)
    m_p_total = mass_paper_g / 1000   # Total paper mass (kg)
    m_c = mass_cylinder_g / 1000      # Cylinder mass (kg)
    
    # --- Derived Initial Properties ---
    # Outer radius of the full roll
    r_o = r_i + num_wraps * tau
    # Total length of the paper, calculated from its cross-sectional area and thickness
    L = np.pi * (r_o**2 - r_i**2) / tau

    print("--- Problem Parameters ---")
    print(f"Gravity (g): {g} m/s^2")
    print(f"Inner Radius (r_i): {r_i:.3f} m")
    print(f"Paper Thickness (tau): {tau:.4f} m")
    print(f"Paper Mass (m_p_total): {m_p_total:.3f} kg")
    print(f"Cylinder Mass (m_c): {m_c:.3f} kg")
    print(f"Calculated Outer Radius (r_o): {r_o:.3f} m")
    print(f"Calculated Total Paper Length (L): {L:.2f} m")
    print("--------------------------\n")
    
    def get_derivatives(state):
        """
        Calculates the derivatives dy/dt and dv/dt for the RK4 solver.
        state[0] = y (distance fallen)
        state[1] = v (velocity)
        """
        y, v = state

        # Cap y at a value just below L to prevent math errors from overshooting.
        if y >= L:
            y = L * (1 - 1e-9)

        # Calculate instantaneous radius, r(y)
        r_sq = r_o**2 - y * tau / np.pi
        r = np.sqrt(r_sq)

        # Calculate instantaneous mass, M(y)
        m_p_remaining = m_p_total * (1.0 - y / L)
        M = m_c + m_p_remaining

        # Calculate instantaneous moment of inertia, I(y)
        I_c = m_c * r_i**2
        I_p = 0.5 * m_p_remaining * (r**2 + r_i**2)
        I = I_c + I_p
        
        # This term is needed for the full, correct equation of motion
        dr_dy = -tau / (2 * np.pi * r)

        # Calculate acceleration a = dv/dt using the full equation of motion
        numerator = g + (I * v**2 / (M * r**3)) * dr_dy
        denominator = 1 + I / (M * r**2)
        a = numerator / denominator

        return np.array([v, a])

    # --- RK4 Solver Initialization ---
    state = np.array([0.0, 0.0])  # Initial state: [y_0, v_0]
    t = 0.0
    dt = 0.001  # Time step in seconds

    # --- RK4 Simulation Loop ---
    while state[0] < L:
        state_prev = state.copy()
        t_prev = t

        k1 = dt * get_derivatives(state)
        k2 = dt * get_derivatives(state + k1/2)
        k3 = dt * get_derivatives(state + k2/2)
        k4 = dt * get_derivatives(state + k3)
        
        state += (k1 + 2*k2 + 2*k3 + k4) / 6
        t += dt

    # --- Final Time Interpolation ---
    # The loop stops when y > L. We interpolate to find the exact time y = L.
    # Fraction of the last time step needed to reach L
    fraction = (L - state_prev[0]) / (state[0] - state_prev[0])
    final_time = t_prev + dt * fraction
    
    # Print the equation used for the final time interpolation for clarity
    print(f"Final time (t) is found using linear interpolation:")
    print(f"t = t_prev + dt * (L - y_prev) / (y_curr - y_prev)")
    print(f"t = {t_prev:.3f} + {dt:.3f} * ({L:.3f} - {state_prev[0]:.3f}) / ({state[0]:.3f} - {state_prev[0]:.3f})")
    print(f"t = {final_time:.2f} seconds\n")
    
    print("Final Answer:")
    print(f"The time it takes for the toilet paper to completely unroll is {final_time:.2f} seconds.")
    return final_time

# Run the simulation and print the final answer
final_time = solve_toilet_paper_fall()
# <<<final_answer>>>
<<<3.17>>>