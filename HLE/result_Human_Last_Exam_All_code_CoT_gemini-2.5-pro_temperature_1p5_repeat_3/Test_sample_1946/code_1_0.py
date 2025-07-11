import numpy as np

def solve_toilet_paper_fall():
    """
    Solves the falling toilet paper problem using a Runge-Kutta numerical method.
    """
    # --- 1. Define Constants and Initial Parameters ---
    g = 9.81  # Acceleration due to gravity (m/s^2)
    m_paper_total = 0.2  # Total mass of the paper (kg)
    m_tube = 0.02  # Mass of the cardboard tube (kg)
    d_tube = 0.04  # Diameter of the tube (m)
    r_tube = d_tube / 2  # Radius of the tube (m)
    thickness = 0.0005  # Thickness of the paper (m)
    N = 100  # Number of wraps
    h = 50.0 # Initial height (m)

    # --- 2. Calculate Derived Properties of the Roll ---
    # Initial outer radius
    r_outer_initial = r_tube + N * thickness
    
    # Total cross-sectional area of the paper
    paper_area = np.pi * (r_outer_initial**2 - r_tube**2)
    
    # Total length of the paper sheet
    L_total = paper_area / thickness

    # --- 3. Define the Function for Instantaneous Acceleration ---
    def get_acceleration(y):
        """
        Calculates the acceleration of the roll given the unrolled length y.
        The acceleration is based on the equation a = g / (1 + I / (m * r^2)).
        """
        # Ensure y does not exceed the total length
        if y >= L_total:
            y = L_total

        # Fraction of paper that has been unrolled and that remains
        unrolled_fraction = y / L_total
        remaining_fraction = 1.0 - unrolled_fraction

        # Calculate the mass, radius, and moment of inertia of the remaining roll
        m_paper_rem = m_paper_total * remaining_fraction
        m_total = m_tube + m_paper_rem

        # If mass is zero (or less), acceleration is zero.
        if m_total <= 1e-9: return 0.0

        rem_area = paper_area * remaining_fraction
        # From rem_area = pi*(r_outer^2 - r_tube^2), we solve for r_outer^2
        r_outer_sq = r_tube**2 + rem_area / np.pi
        
        # Moment of inertia for the hollow cylinder of paper
        I_paper = 0.5 * m_paper_rem * (r_outer_sq + r_tube**2)
        # Moment of inertia for the cardboard tube
        I_tube = m_tube * r_tube**2
        # Total moment of inertia
        I_total = I_tube + I_paper

        # The final acceleration formula
        a = g / (1 + I_total / (m_total * r_outer_sq))
        
        return a

    # --- 4. RK4 Numerical Solver ---
    t = 0.0      # Initial time
    y = 0.0      # Initial position (unrolled length)
    v = 0.0      # Initial velocity
    dt = 0.001   # Time step in seconds

    # Store previous state for final interpolation
    t_prev, y_prev = t, y

    # Run the simulation until the paper is fully unrolled
    while y < L_total:
        t_prev, y_prev = t, y

        # RK4 steps for the ODE system {dy/dt = v, dv/dt = a(y)}
        k1_v = v
        k1_a = get_acceleration(y)
        
        k2_v = v + dt / 2.0 * k1_a
        k2_a = get_acceleration(y + dt / 2.0 * k1_v)

        k3_v = v + dt / 2.0 * k2_a
        k3_a = get_acceleration(y + dt / 2.0 * k2_v)
        
        k4_v = v + dt * k3_a
        k4_a = get_acceleration(y + dt * k3_v)

        # Update position, velocity, and time
        y += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        v += (dt / 6.0) * (k1_a + 2*k2_a + 2*k3_a + k4_a)
        t += dt

    # --- 5. Interpolate for Final Time and Print Results ---
    # The loop stops when y > L_total. Linearly interpolate between the last two points 
    # (t_prev, y_prev) and (t, y) to find the exact time when y = L_total.
    time_to_unroll = t_prev + (L_total - y_prev) * (t - t_prev) / (y - y_prev)

    # Print a summary of the problem and the governing equation
    print("--- Toilet Paper Unrolling Problem ---")
    print(f"Total length of paper: {L_total:.2f} m (which is less than the initial height of {h:.0f} m)")
    print("The toilet paper will completely unroll before hitting the ground.")
    
    # Calculate initial values for the equation showcase
    m0 = m_tube + m_paper_total
    r0 = r_outer_initial
    I0 = (m_tube * r_tube**2) + (0.5 * m_paper_total * (r0**2 + r_tube**2))
    a0 = get_acceleration(0)

    print("\nThe governing equation for acceleration (a) is: a = g / (1 + I / (m * r^2))")
    print("where g is gravity, I is moment of inertia, m is mass, and r is the outer radius.")
    print("These values change as the roll unrolls. For example, at the start (t=0):")
    print(f"  g = {g} m/s^2")
    print(f"  I = {I0:.6f} kg*m^2")
    print(f"  m = {m0:.3f} kg")
    print(f"  r = {r0:.3f} m")
    print(f"This gives an initial acceleration a(0) = {a0:.2f} m/s^2.\n")

    print(f"The simulation shows it takes the following time for the roll to completely unravel:")
    print(f"{time_to_unroll:.2f} seconds")

    return time_to_unroll

# Execute the solver
final_time = solve_toilet_paper_fall()
print(f"<<<{final_time:.2f}>>>")
