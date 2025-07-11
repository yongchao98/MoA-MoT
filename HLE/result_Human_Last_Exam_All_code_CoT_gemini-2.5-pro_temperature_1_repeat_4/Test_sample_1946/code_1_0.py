import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a falling toilet paper roll to completely unroll.

    The calculation is done by numerically solving the equations of motion using the
    4th-order Runge-Kutta (RK4) method.
    """

    # 1. Define physical constants from the problem description
    g = 9.81  # Gravitational acceleration, m/s^2
    m_c = 0.020  # Mass of the cardboard cylinder, kg
    d_c = 0.04  # Diameter of the cardboard cylinder, m
    m_p_total = 0.200  # Total mass of the paper, kg
    paper_thickness = 0.0005  # Thickness of the paper, m (0.5 mm)
    paper_wraps = 100  # Number of times the paper is wrapped around the roll

    # 2. Calculate derived constants
    r_c = d_c / 2  # Radius of the cardboard cylinder, m
    # Initial outer radius of the full roll
    r_0 = r_c + paper_wraps * paper_thickness
    # Total length of the paper, calculated as the sum of circumferences of all wraps
    # Using the formula for the sum of an arithmetic series: N * pi * (r_first + r_last)
    # This simplifies to N * pi * (2*r_c + (N+1)*t) if we average layer radii
    # A common approximation is L = N * pi * (2*r_c + N*t)
    L = paper_wraps * np.pi * (2 * r_c + paper_wraps * paper_thickness)

    # Print the parameters used in the model
    print("--- Model Parameters ---")
    print(f"Gravitational acceleration (g): {g} m/s^2")
    print(f"Cardboard cylinder mass (m_c): {m_c} kg")
    print(f"Cardboard cylinder radius (r_c): {r_c} m")
    print(f"Total paper mass (M_p): {m_p_total} kg")
    print(f"Paper thickness (t): {paper_thickness} m")
    print(f"Number of paper wraps (N): {paper_wraps}")
    print(f"Initial roll radius (R_0): {r_0:.4f} m")
    print(f"Total paper length (L): {L:.4f} m")
    print("------------------------")

    # 3. Define functions for state-dependent properties (as a function of unrolled length l)
    def get_radius(l):
        """Calculates the outer radius of the roll when length l has unrolled."""
        if l >= L:
            return r_c
        # The cross-sectional area of the paper is proportional to the remaining length.
        # R(l)^2 = r_c^2 + (r_0^2 - r_c^2) * (remaining_length / total_length)
        radius_sq = r_c**2 + (r_0**2 - r_c**2) * (1 - l / L)
        return np.sqrt(radius_sq)

    def get_mass(l):
        """Calculates the total mass of the roll (cylinder + remaining paper)."""
        if l >= L:
            return m_c
        paper_mass = m_p_total * (1 - l / L)
        return m_c + paper_mass

    def get_moment_of_inertia(l):
        """Calculates the moment of inertia of the roll."""
        # Moment of inertia of the cardboard tube (thin hollow cylinder: I = mr^2)
        I_c = m_c * r_c**2
        if l >= L:
            return I_c
        
        # Moment of inertia of the remaining paper (hollow cylinder: I = 1/2 * m * (r1^2 + r2^2))
        m_p = m_p_total * (1 - l / L)
        R_l = get_radius(l)
        I_p = 0.5 * m_p * (R_l**2 + r_c**2)
        return I_c + I_p

    # 4. Define the derivative function for the ODE system dy/dt = f(t, y)
    # The state vector is y = [l, v], where l is unrolled length and v is velocity.
    def f(t, y):
        l, v = y[0], y[1]

        if l >= L:
            return np.array([0.0, 0.0])

        m_l = get_mass(l)
        R_l = get_radius(l)
        I_l = get_moment_of_inertia(l)
        
        # The acceleration 'a' is derived from F=ma and tau=I*alpha, with the
        # kinematic constraint v = omega * R. Differentiating this constraint
        # introduces a term related to the change in radius.
        # a = (g*m*R^2 - I*v^2*(R_0^2-r_c^2)/(2*L*R^2)) / (m*R^2 + I)
        
        numerator = g * m_l * R_l**2 - I_l * v**2 * (r_0**2 - r_c**2) / (2 * L * R_l**2)
        denominator = m_l * R_l**2 + I_l
        
        a = numerator / denominator if denominator else 0.0
        
        return np.array([v, a])

    # 5. Implement the RK4 solver
    dt = 0.001  # Time step for the simulation, in seconds
    t = 0.0
    # Initial state: l=0 (no paper unrolled), v=0 (starts from rest)
    y = np.array([0.0, 0.0])
    y_prev, t_prev = y, t

    while y[0] < L:
        y_prev, t_prev = y, t
        
        k1 = f(t, y)
        k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1)
        k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2)
        k4 = f(t + dt, y + dt * k3)
        
        y = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t = t + dt

    # 6. Interpolate to find the exact time when l = L
    l_prev = y_prev[0]
    l_curr = y[0]
    
    # t_final = t_prev + dt * (L - l_prev) / (l_curr - l_prev)
    if (l_curr - l_prev) == 0:
        t_final = t
    else:
        t_final = t_prev + (t - t_prev) * (L - l_prev) / (l_curr - l_prev)

    print(f"\nCalculation finished.")
    print(f"The time for the toilet paper to unroll completely is: {t_final:.2f} seconds.")
    return t_final

# Execute the solver and print the final answer in the required format
final_time = solve_toilet_paper_fall()
print(f"\n<<< {final_time:.2f} >>>")
