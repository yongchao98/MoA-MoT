import math

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to unroll completely when dropped.
    The calculation is performed by numerically solving the equation of motion using the RK4 method.
    """

    # 1. Define physical constants from the problem description in SI units.
    g = 9.81  # Gravitational acceleration (m/s^2)
    d_inner = 0.04  # Inner cardboard cylinder diameter (m)
    t_paper = 0.0005 # Toilet paper thickness (m)
    N_wraps = 100   # Number of times paper is wrapped around the roll
    M_p = 0.200     # Mass of paper (kg)
    M_c = 0.020     # Mass of cardboard cylinder (kg)

    # 2. Calculate derived properties of the roll.
    r_inner = d_inner / 2
    r_outer = r_inner + N_wraps * t_paper
    
    # The total length L is the sum of the circumferences of each wrap.
    # L = sum_{n=1 to 100} 2*pi*(r_inner + n*t_paper)
    l_sum_term = N_wraps * r_inner + t_paper * N_wraps * (N_wraps + 1) / 2
    L = 2 * math.pi * l_sum_term
    
    # Moment of inertia of the hollow cardboard cylinder (constant)
    I_c = M_c * r_inner**2

    print("--- System Parameters ---")
    print(f"Paper Mass (M_p):      {M_p:.3f} kg")
    print(f"Core Mass (M_c):       {M_c:.3f} kg")
    print(f"Inner Radius (r_i):    {r_inner:.3f} m")
    print(f"Outer Radius (R_0):    {r_outer:.3f} m")
    print(f"Total Paper Length (L):{L:.2f} m")
    print(f"Gravity (g):           {g:.2f} m/s^2")
    print("--------------------------\n")
    
    # 3. Define functions for properties that change with unrolled length y.
    
    def get_roll_properties(y):
        """Calculates mass, radius squared, and moment of inertia for a given unrolled length y."""
        # Fraction of paper remaining
        frac_remaining = (1 - y / L) if y < L else 0

        # Current mass of the roll
        current_mass = M_c + M_p * frac_remaining
        
        # Current outer radius squared of the paper roll
        current_R_sq = r_inner**2 + (r_outer**2 - r_inner**2) * frac_remaining
        
        # Current moment of inertia of the paper part
        m_p_y = M_p * frac_remaining
        I_p_y = 0.5 * m_p_y * (current_R_sq + r_inner**2)
        current_I = I_c + I_p_y
        
        return current_mass, current_R_sq, current_I

    def get_acceleration(y):
        """Calculates acceleration a(y) based on the unrolled length y."""
        if y >= L:
            return g / 2.0  # Theoretical limit as the roll becomes just the core
            
        m_y, R_sq_y, I_y = get_roll_properties(y)
        
        if m_y <= 0 or R_sq_y <= 0: return 0

        # This factor is k from the yo-yo acceleration equation a = g / (1 + k)
        # where k = I / (m*R^2)
        factor = I_y / (m_y * R_sq_y)
        
        return g / (1 + factor)

    # 4. Solve the ODE using the 4th-order Runge-Kutta (RK4) method.
    y, v, t = 0.0, 0.0, 0.0  # Initial state: y(unrolled length), v(velocity), t(time)
    dt = 0.001  # Time step in seconds

    y_prev, t_prev = y, t

    while y < L:
        y_prev, t_prev = y, t

        # RK4 step
        k1_v = dt * get_acceleration(y)
        k1_y = dt * v
        
        k2_v = dt * get_acceleration(y + 0.5 * k1_y)
        k2_y = dt * (v + 0.5 * k1_v)

        k3_v = dt * get_acceleration(y + 0.5 * k2_y)
        k3_y = dt * (v + 0.5 * k2_v)

        k4_v = dt * get_acceleration(y + k3_y)
        k4_y = dt * (v + k3_v)

        # Update state variables
        y += (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6.0
        v += (k1_v + 2*k2_v + 2*k3_v + k4_v) / 6.0
        t += dt

    # 5. Interpolate for the precise time when y reached L.
    # The loop stops when y > L, so we correct for the overshoot in the last step.
    overshoot_dist = y - L
    last_step_dist = y - y_prev
    time_correction = dt * (last_step_dist - overshoot_dist) / last_step_dist
    final_time = t_prev + time_correction

    print(f"Time to unroll completely: {final_time:.2f} seconds.")

# Run the simulation
solve_toilet_paper_fall()