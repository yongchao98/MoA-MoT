import math

def solve_toilet_paper_problem():
    """
    This function solves the falling toilet paper problem by numerically integrating
    the equations of motion using the RK4 method.
    """
    # --- Step 1: Define physical constants from the problem ---
    # Gravitational acceleration (m/s^2)
    G = 9.81
    # Inner cylinder diameter is 4 cm, so radius is 2 cm (m)
    R_INNER = 0.02
    # Paper thickness is 0.5 mm (m)
    PAPER_THICKNESS = 0.0005
    # Number of wraps
    N_WRAPS = 100
    # Mass of the paper (kg)
    M_PAPER = 0.2
    # Mass of the cardboard cylinder (kg)
    M_CYLINDER = 0.02

    # --- Step 2: Derive expressions for changing properties ---
    # Outer radius of the full roll (m)
    R_OUTER = R_INNER + N_WRAPS * PAPER_THICKNESS
    # Total length of the paper, calculated using the conserved volume/area method.
    # This is more accurate than summing circumferences for a continuous sheet.
    # Area = pi * (R_outer^2 - R_inner^2) = L * PAPER_THICKNESS
    L = math.pi * (R_OUTER**2 - R_INNER**2) / PAPER_THICKNESS

    # --- Step 3: Define the acceleration function a(l) ---
    def get_acceleration(l):
        """Calculates the acceleration of the roll for a given unrolled length 'l'."""
        if l >= L:
            return 0

        # Fraction of paper area/mass/volume remaining
        frac_rem = 1.0 - l / L
        
        # Instantaneous mass of remaining paper
        m_p_rem = M_PAPER * frac_rem
        # Instantaneous total mass of the falling roll
        m_total = M_CYLINDER + m_p_rem
        
        # Instantaneous outer radius squared of the paper roll
        r_sq = R_INNER**2 + frac_rem * (R_OUTER**2 - R_INNER**2)
        
        # Instantaneous moment of inertia of the roll
        # I_cylinder = M_CYLINDER * R_INNER^2 (approximating as a thin ring)
        i_cyl = M_CYLINDER * R_INNER**2
        # I_paper = 1/2 * m_p_rem * (r(l)^2 + R_INNER^2) (for a thick hollow cylinder)
        i_paper = 0.5 * m_p_rem * (r_sq + R_INNER**2)
        i_total = i_cyl + i_paper

        # Avoid division by zero if radius becomes extremely small
        if r_sq < 1e-9:
            return 0
        
        # The acceleration is given by a = M*g / (M + I/r^2)
        acceleration = (m_total * G) / (m_total + i_total / r_sq)
        return acceleration

    # --- Step 4: Solve the ODE using RK4 numerical method ---
    t = 0.0       # time
    l = 0.0       # unrolled length
    v = 0.0       # velocity
    dt = 0.001    # time step in seconds, for accuracy

    l_prev, t_prev = 0, 0

    while l < L:
        l_prev, t_prev = l, t
        
        # Standard RK4 implementation for the system:
        # y' = f(y), where y = [l, v]
        k1_l, k1_v = v, get_acceleration(l)
        k2_l, k2_v = v + 0.5*dt*k1_v, get_acceleration(l + 0.5*dt*k1_l)
        k3_l, k3_v = v + 0.5*dt*k2_v, get_acceleration(l + 0.5*dt*k2_l)
        k4_l, k4_v = v + dt*k3_v, get_acceleration(l + dt*k3_l)
        
        # Update state variables l (length) and v (velocity)
        l += (dt/6.0) * (k1_l + 2*k2_l + 2*k3_l + k4_l)
        v += (dt/6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        t += dt

    # The loop stops when l > L. Interpolate for a more accurate final time.
    if l > l_prev:
        t_final = t_prev + dt * (L - l_prev) / (l - l_prev)
    else:
        t_final = t

    # --- Final Output ---
    # The prompt asks to "output each number in the final equation".
    # The final result is not from a single simple equation, but the result of the simulation.
    # We will state the governing equation and the values of the constants used.
    print("The time 't' is found by solving the differential equation d²l/dt² = a(l).")
    print("The acceleration a(l) is given by the equation: a(l) = M(l)*g / (M(l) + I(l)/r(l)²)\n")
    print("The constant values used in this equation are:")
    print(f"g (Gravitational Acceleration) = {G} m/s²")
    print(f"L (Total Paper Length) = {L:.4f} m")
    print(f"m_p (Initial Paper Mass) = {M_PAPER} kg")
    print(f"m_c (Cylinder Mass) = {M_CYLINDER} kg")
    print(f"r_out (Initial Radius) = {R_OUTER} m")
    print(f"r_in (Final Radius) = {R_INNER} m\n")

    print("The final calculated time to unroll completely is:")
    print(f"t = {t_final:.2f} seconds")
    
    return t_final

if __name__ == '__main__':
    final_time = solve_toilet_paper_problem()
    # The final answer is returned in the specified format.
    print(f"\n<<<{final_time:.2f}>>>")