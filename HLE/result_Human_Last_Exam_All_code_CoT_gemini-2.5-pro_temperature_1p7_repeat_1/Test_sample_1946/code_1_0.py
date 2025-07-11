import math

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to unroll completely
    when dropped, using the principle of energy conservation and numerical integration.
    """

    # 1. Define constants from the problem statement
    g = 9.8  # m/s^2
    d_inner = 4.0  # cm
    r_inner = d_inner / 2 / 100  # convert cm to m
    paper_thickness = 0.5 / 1000  # convert mm to m
    N_wraps = 100
    M_paper = 200.0 / 1000  # convert g to kg
    M_cardboard = 20.0 / 1000  # convert g to kg

    # 2. Calculate derived properties of the full roll
    r_outer = r_inner + N_wraps * paper_thickness
    
    # Total length of the paper, calculated from the area of the paper annulus
    # Area = pi * (r_outer^2 - r_inner^2) = L_paper * paper_thickness
    if paper_thickness > 0:
        L_paper = math.pi * (r_outer**2 - r_inner**2) / paper_thickness
    else:
        L_paper = 0

    if L_paper == 0:
        print("Paper has no length. Time is 0.")
        return

    # Linear mass density of the paper
    lambda_paper = M_paper / L_paper

    # Moment of inertia of the cardboard tube (approximated as a thin cylinder)
    I_cardboard = M_cardboard * r_inner**2

    # 3. Define a function to calculate the squared velocity v^2(x)
    #    based on the unrolled length x, using energy conservation (KE + PE = 0).
    def get_v_squared(x):
        # PE = -g * (m_roll * x + m_unrolled * x / 2)
        # KE = v^2 * (0.5 * m_roll + 0.5 * I_x / r_sq_x + (1/6) * m_unrolled)
        # From PE + KE = 0, we solve for v^2.

        if x <= 1e-9:  # Avoid division by zero at the start
            return 0

        # Properties as a function of unrolled length x
        m_unrolled = lambda_paper * x
        m_paper_rem = M_paper - m_unrolled
        m_roll = M_cardboard + m_paper_rem

        # Current outer radius squared of the roll
        r_sq_x = r_inner**2 + (L_paper - x) * paper_thickness / math.pi
        
        # Current moment of inertia of the roll
        I_paper_rem = 0.5 * m_paper_rem * (r_sq_x + r_inner**2)
        I_x = I_cardboard + I_paper_rem

        # Numerator of v^2 from Potential Energy release
        pe_release = g * (m_roll * x + m_unrolled * x / 2.0)
        
        # Denominator of v^2 from Kinetic Energy terms
        ke_factor = 0.5 * m_roll + 0.5 * I_x / r_sq_x + (1.0 / 6.0) * m_unrolled

        if ke_factor == 0:
            return 0
        
        return pe_release / ke_factor

    # 4. Numerically integrate dt = dx / v(x) from x=0 to x=L
    N_steps = 2000000  # Number of steps for the integration (high for accuracy)
    dx = L_paper / N_steps
    total_time = 0.0

    # We integrate using the midpoint rule for better accuracy.
    # The singularity at x=0 (where v=0) is integrable.
    for i in range(1, N_steps + 1):
        x_i = (i - 0.5) * dx
        v_sq = get_v_squared(x_i)
        
        if v_sq > 0:
            v_i = math.sqrt(v_sq)
            dt = dx / v_i
            total_time += dt

    # 5. Output the parameters and the final answer
    print("This script calculates the time for a toilet paper roll to unravel when dropped.")
    print("\nParameters used in the calculation:")
    print(f"  Inner cylinder diameter: {d_inner} cm")
    print(f"  Paper thickness: {paper_thickness*1000} mm")
    print(f"  Number of wraps: {N_wraps}")
    print(f"  Mass of paper: {M_paper*1000} g")
    print(f"  Mass of cardboard cylinder: {M_cardboard*1000} g")
    print(f"  Gravitational acceleration: {g} m/s^2")

    print("\nDerived value:")
    print(f"  Total length of paper: {L_paper:.2f} m")

    print("\nResult:")
    print(f"Time for the toilet paper to completely unroll: {total_time:.2f} s")

# Run the simulation
solve_toilet_paper_fall()