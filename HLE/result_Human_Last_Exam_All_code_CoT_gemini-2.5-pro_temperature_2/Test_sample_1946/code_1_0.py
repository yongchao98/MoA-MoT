import math

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to completely unroll while falling.
    This is solved by numerically integrating the equations of motion using the RK4 method.
    """
    # 1. Define constants and given parameters
    g = 9.81  # Gravitational acceleration in m/s^2

    # Physical properties of the toilet paper roll
    d_inner_cm = 4.0  # cm
    t_paper_mm = 0.5  # mm
    n_wraps = 100
    m_paper_total_g = 200  # grams
    m_tube_g = 20  # grams

    # Convert to SI units (meters, kilograms)
    r_inner = (d_inner_cm / 100) / 2.0
    t_paper = t_paper_mm / 1000
    m_paper_total = m_paper_total_g / 1000
    m_tube = m_tube_g / 1000

    # 2. Calculate derived parameters
    # Initial outer radius of the full roll
    r_outer_initial = r_inner + n_wraps * t_paper

    # Total length of the paper, L, derived from the area
    # Area = pi * (R_outer^2 - R_inner^2) = L * t
    paper_area = math.pi * (r_outer_initial**2 - r_inner**2)
    L = paper_area / t_paper

    # Moment of inertia of the cardboard tube, modeled as a thin hoop (I = mr^2)
    I_tube = m_tube * r_inner**2

    # 3. Define the function for the system's acceleration
    def get_acceleration(y):
        """Calculates acceleration 'a' for a given unrolled length 'y'."""
        if y >= L:
            return 0  # Stop calculation when fully unrolled

        # Calculate remaining paper fraction and mass
        rem_paper_fraction = (L - y) / L
        m_falling = m_tube + m_paper_total * rem_paper_fraction
        
        # The roll has unrolled completely; it's now just the tube.
        # This part of the code might be needed if dt is very large, but the loop will terminate anyway.
        if rem_paper_fraction <= 1e-9:
             r_outer_sq = r_inner**2
             I_paper = 0
        else:
            # Calculate current outer radius squared of the roll
            r_outer_sq = r_inner**2 + (L - y) * t_paper / math.pi
            
            # Calculate moment of inertia of the remaining paper (thick hollow cylinder)
            m_paper_rem = m_paper_total * rem_paper_fraction
            I_paper = 0.5 * m_paper_rem * (r_outer_sq + r_inner**2)
            
        # Total moment of inertia
        I_total = I_tube + I_paper

        # Calculate acceleration using the formula a = Mg / (M + I/R^2)
        acceleration = (m_falling * g) / (m_falling + I_total / r_outer_sq)
        
        return acceleration

    # 4. Set up and run the RK4 numerical solver
    y = 0.0  # Initial unrolled length
    v = 0.0  # Initial velocity
    t = 0.0  # Initial time
    dt = 0.001  # Time step in seconds

    y_prev, t_prev = y, t

    # Loop until the paper is fully unrolled (y >= L)
    while y < L:
        t_prev, y_prev = t, y

        # Standard RK4 algorithm for the system dy/dt = v and dv/dt = a(y)
        k1_v = v
        k1_a = get_acceleration(y)
        
        k2_v = v + 0.5 * dt * k1_a
        k2_a = get_acceleration(y + 0.5 * dt * k1_v)

        k3_v = v + 0.5 * dt * k2_a
        k3_a = get_acceleration(y + 0.5 * dt * k2_v)

        k4_v = v + dt * k3_a
        k4_a = get_acceleration(y + dt * k3_v)

        # Update position (y) and velocity (v)
        y += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        v += (dt / 6.0) * (k1_a + 2*k2_a + 2*k3_a + k4_a)
        
        # Update time
        t += dt

    # 5. Interpolate for a more precise final time
    # When the loop finishes, y is slightly > L. Find the time t when y was exactly L.
    # final_time = t_prev + (time_to_cover_remaining_dist)
    # t_final = t_prev + (L - y_prev) / (v_avg) ~= t_prev + dt * (L - y_prev) / (y - y_prev)
    final_time = t_prev + dt * (L - y_prev) / (y - y_prev)
    
    print(f"Total length of the paper is: {L:.2f} m.")
    print("The final time is found by interpolating the last step of the simulation:")
    print(f"Final time = {t_prev:.4f}s + ({dt:.3f}s * ({L:.4f}m - {y_prev:.4f}m) / ({y:.4f}m - {y_prev:.4f}m))")
    print(f"The time it takes for the toilet paper to completely unroll is {final_time:.2f} seconds.")
    print(f"<<<{final_time:.2f}>>>")

solve_toilet_paper_fall()