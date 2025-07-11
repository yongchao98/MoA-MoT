import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a toilet paper roll to fall and unroll completely.
    
    This function models the physics of the falling roll using Lagrangian mechanics
    and solves the resulting differential equation of motion using the 4th-order
    Runge-Kutta method.
    """
    # 1. Define the physical parameters given in the problem
    g = 9.8          # Acceleration due to gravity in m/s^2
    d_c = 0.04       # Diameter of the cardboard cylinder in m
    m_c = 0.02       # Mass of the cardboard cylinder in kg
    t_paper = 0.0005 # Thickness of a single sheet of paper in m
    N_wraps = 100    # Number of times the paper is wrapped around the roll
    M_p = 0.2        # Total mass of the paper in kg
    
    # 2. Calculate derived parameters based on the initial values
    r_c = d_c / 2.0  # Radius of the cardboard cylinder

    # We use a continuous model to find the total length of the paper.
    # The initial outer radius of the full roll:
    R_outer = r_c + N_wraps * t_paper
    # The total length of the unrolled paper:
    L = np.pi * (R_outer**2 - r_c**2) / t_paper
    # The linear mass density of the paper (mass per unit length):
    lambda_p = M_p / L

    # This constant represents how quickly the roll's radius squared changes as it unrolls.
    # It is the derivative d(R^2)/dy.
    C_R_sq = -(R_outer**2 - r_c**2) / L

    # 3. Define the system of differential equations derived from the Lagrangian
    def get_acceleration(y, v):
        """
        Calculates the acceleration of the roll at a given position y and velocity v.
        y: the length of paper that has unrolled (m).
        v: the instantaneous velocity of the roll (m/s).
        """
        # The radius squared of the paper roll as a function of unrolled length y
        R_sq = r_c**2 + (R_outer**2 - r_c**2) / L * (L - y)

        # The effective mass M_eff(y) includes both translational and rotational inertia
        m_eff_trans = m_c * (1 + r_c**2 / R_sq)
        m_eff_rot = lambda_p * (L - y) * (1.5 + 0.5 * r_c**2 / R_sq)
        M_eff = m_eff_trans + m_eff_rot

        # The derivative of the effective mass with respect to y, M_eff'(y)
        d_inv_R_sq_dy = -C_R_sq / (R_sq**2)
        
        m_prime_term1 = m_c * r_c**2 * d_inv_R_sq_dy
        
        # Product rule for the second term of M_eff
        prod1 = -lambda_p * (1.5 + 0.5 * r_c**2 / R_sq)
        prod2 = lambda_p * (L - y) * (0.5 * r_c**2 * d_inv_R_sq_dy)
        m_prime_term2 = prod1 + prod2
        
        M_prime = m_prime_term1 + m_prime_term2

        # The effective gravitational force term from the Lagrangian equation
        F_grav_eff = g * (m_c + lambda_p * (L - y))

        # The final equation for acceleration a = dv/dt
        numerator = F_grav_eff - 0.5 * v**2 * M_prime
        a = numerator / M_eff
        
        return a

    # 4. Implement the RK4 numerical solver
    y = 0.0          # Initial unrolled length
    v = 0.0          # Initial velocity
    t = 0.0          # Initial time
    dt = 0.001       # Time step for the simulation in seconds

    y_prev, t_prev = 0.0, 0.0

    while y < L:
        y_prev, t_prev = y, t

        # Standard RK4 algorithm
        k1_v = get_acceleration(y, v)
        k1_y = v
        
        k2_v = get_acceleration(y + 0.5 * dt * k1_y, v + 0.5 * dt * k1_v)
        k2_y = v + 0.5 * dt * k1_v
        
        k3_v = get_acceleration(y + 0.5 * dt * k2_y, v + 0.5 * dt * k2_v)
        k3_y = v + 0.5 * dt * k2_v
        
        k4_v = get_acceleration(y + dt * k3_y, v + dt * k3_v)
        k4_y = v + dt * k3_v

        # Update velocity and position
        v += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        y += (dt / 6.0) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
        t += dt

    # 5. Interpolate for a precise final time when y exactly equals L
    if y > y_prev:
        t_final = t_prev + dt * (L - y_prev) / (y - y_prev)
    else:
        t_final = t

    # 6. Print the results and explanation
    print("This solution uses the following physical parameters:")
    print(f"  - Cardboard cylinder diameter: {d_c*100:.0f} cm")
    print(f"  - Cardboard cylinder mass: {m_c*1000:.0f} g")
    print(f"  - Paper thickness: {t_paper*1000} mm")
    print(f"  - Number of wraps: {N_wraps}")
    print(f"  - Total paper mass: {M_p*1000:.0f} g")
    print("\nDerived parameters for the model:")
    print(f"  - Total paper length (L): {L:.2f} m")
    print(f"  - Paper linear mass density (λ): {lambda_p:.5f} kg/m")
    
    print("\nThe equation for the roll's acceleration is solved numerically.")
    print("\nResult:")
    print(f"The time it takes for the toilet paper to completely unroll is {t_final:.2f} seconds.")
    
    # Final answer in the required format
    print(f"\n<<<{t_final:.2f}>>>")

# Run the simulation
solve_toilet_paper_fall()