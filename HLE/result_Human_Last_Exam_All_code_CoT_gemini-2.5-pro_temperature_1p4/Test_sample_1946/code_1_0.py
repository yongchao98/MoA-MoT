import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time for a falling toilet paper roll to completely unroll.

    This function models the physics using the work-energy theorem and solves the
    resulting differential equations using the 4th-order Runge-Kutta method.
    """
    # 1. Physical Parameters from the problem statement
    g = 9.8  # Acceleration due to gravity (m/s^2)
    d_cyl = 4.0 / 100  # Cardboard cylinder diameter (m)
    r_cyl = d_cyl / 2.0  # Cardboard cylinder radius (m)
    t_paper = 0.5 / 1000  # Paper thickness (m)
    n_wraps = 100  # Number of wraps
    m_paper = 200.0 / 1000  # Total mass of paper (kg)
    m_cyl = 20.0 / 1000  # Mass of cardboard cylinder (kg)

    # 2. Derived Geometric and Mass Properties
    r_outer = r_cyl + n_wraps * t_paper
    # Total paper length calculated from cross-sectional area
    l_total = np.pi * (r_outer**2 - r_cyl**2) / t_paper
    # Linear mass density of the paper
    lam = m_paper / l_total
    # A constant needed for derivatives
    dS_dy_const = -(r_outer**2 - r_cyl**2) / l_total

    print("--- System Parameters ---")
    print(f"Cardboard mass (Mc): {m_cyl} kg")
    print(f"Paper mass (Mp): {m_paper} kg")
    print(f"Inner radius (rcyl): {r_cyl} m")
    print(f"Outer radius (Router): {r_outer:.3f} m")
    print(f"Total paper length (L_total): {l_total:.2f} m")
    print("-" * 25)

    def get_acceleration(y, v):
        """
        Calculates acceleration a(y, v) from the energy conservation equation.
        This function is the core of the ODE system.
        """
        # Ensure y does not exceed the total length for calculations
        y_calc = min(y, l_total)

        # Instantaneous radius squared (S)
        r_sq = r_cyl**2 + (r_outer**2 - r_cyl**2) * (1 - y_calc / l_total)
        
        # Instantaneous mass of rolled paper
        mp_rolled = m_paper * (1 - y_calc / l_total)
        
        # Effective mass term from the kinetic energy equation (translational + rotational)
        # M_eff = m_falling + I_falling / r^2
        m_eff = m_cyl * (1 + r_cyl**2 / (2 * r_sq)) + mp_rolled * (1.5 + r_cyl**2 / (2 * r_sq))

        if m_eff <= 0: return 0

        # Derivative of the effective mass w.r.t y (for the acceleration equation)
        dMeff_dy_term1 = (-r_cyl**2 / (2 * r_sq**2)) * dS_dy_const * (m_cyl + mp_rolled)
        dMeff_dy_term2 = -lam * (1.5 + r_cyl**2 / (2 * r_sq))
        dMeff_dy = dMeff_dy_term1 + dMeff_dy_term2
        
        # Acceleration a = [ d(PE)/dy - 0.5*v^2*d(M_eff)/dy ] / M_eff
        # where d(PE)/dy is the net gravitational force term.
        numerator_a = g * (m_cyl + m_paper - 3 * lam * y_calc) - 0.5 * v**2 * dMeff_dy
        a = numerator_a / m_eff
        
        return a

    # 3. Numerical Integration (RK4)
    t = 0.0
    y = 0.0  # Initial distance fallen
    v = 0.0  # Initial velocity
    dt = 0.001  # Time step for simulation (s)

    a_initial = get_acceleration(0, 0)
    print(f"The equation for initial acceleration is a = g * (Mc + Mp) / M_eff(0)")
    print(f"Initial acceleration: a(0,0) = {a_initial:.2f} m/s^2")
    print("-" * 25)


    y_prev, t_prev = y, t
    while y < l_total:
        y_prev, t_prev = y, t

        # RK4 calculation steps
        k1_v = get_acceleration(y, v)
        k1_y = v

        k2_v = get_acceleration(y + 0.5*dt*k1_y, v + 0.5*dt*k1_v)
        k2_y = v + 0.5*dt*k1_v

        k3_v = get_acceleration(y + 0.5*dt*k2_y, v + 0.5*dt*k2_v)
        k3_y = v + 0.5*dt*k2_v
        
        k4_v = get_acceleration(y + dt*k3_y, v + dt*k3_v)
        k4_y = v + dt*k3_v
        
        # Update velocity and position
        v += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        y += (dt / 6.0) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
        t += dt

    # 4. Interpolate for the final time to get a precise result
    # Find the time when y exactly equals l_total
    if y > y_prev:
        final_time = t_prev + dt * (l_total - y_prev) / (y - y_prev)
    else:
        final_time = t
    
    print("\n--- Final Result ---")
    print(f"The final equation to solve is finding T such that y(T) = L_total.")
    print(f"The time it takes for the toilet paper to fully unroll is:")
    print(f"{final_time:.2f} seconds")
    return final_time

if __name__ == '__main__':
    final_time = solve_toilet_paper_fall()
    # The final answer in the required format
    # print(f"\n<<<{final_time:.2f}>>>")