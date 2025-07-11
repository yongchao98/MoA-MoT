import math
import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to completely
    unroll when dropped from a height, using a numerical ODE solver.
    """
    # Step 1: Define physical constants from the problem statement
    g = 9.8  # Acceleration due to gravity (m/s^2)
    m_paper_total = 0.200  # Total mass of paper (kg)
    m_cylinder = 0.020  # Mass of cardboard cylinder (kg)
    d_cylinder = 0.04  # Diameter of cylinder (m)
    d_paper = 0.0005  # Thickness of paper (m, 0.5 mm)
    n_wraps = 100  # Number of wraps

    # Step 2: Calculate derived parameters of the roll
    r_inner = d_cylinder / 2.0
    r_outer = r_inner + n_wraps * d_paper
    
    r_inner_sq = r_inner**2
    r_outer_sq = r_outer**2

    # Total length of the paper (L) calculated from its volume
    paper_area = math.pi * (r_outer_sq - r_inner_sq)
    L = paper_area / d_paper

    # Linear mass density of the paper (mass per unit length)
    lambda_paper = m_paper_total / L

    # Moment of inertia of the cylinder (modeled as a thin shell)
    I_cylinder = m_cylinder * r_inner_sq
    
    # A constant for calculating remaining paper mass
    C1 = m_paper_total / (r_outer_sq - r_inner_sq)

    # Step 3: Define helper functions for properties that change with distance 'y'
    def get_r_sq(y):
        r_sq = r_outer_sq - y * d_paper / math.pi
        return max(r_inner_sq, r_sq)

    def get_mass_paper_rem(r_sq):
        return C1 * (r_sq - r_inner_sq)

    def get_total_mass_roll(r_sq):
        return m_cylinder + get_mass_paper_rem(r_sq)

    def get_moment_of_inertia_roll(r_sq, m_paper_rem):
        I_paper_rem = 0.5 * m_paper_rem * (r_sq + r_inner_sq)
        return I_cylinder + I_paper_rem

    def get_dK_dy(r_sq):
        if r_sq <= r_inner_sq:
            return 0.0
        dr_sq_dy = -d_paper / math.pi
        r_sq_sq = r_sq**2
        dK_dr_sq = (1.5 * C1) - (I_cylinder / r_sq_sq) + (0.5 * C1 * (r_inner_sq**2) / r_sq_sq)
        return dK_dr_sq * dr_sq_dy

    # Step 4: Define the function for the system of ODEs
    # State S = [y, v], returns dS/dt = [v, a]
    def odefunc(t, S):
        y, v = S[0], S[1]
        if y >= L:
            return np.array([0.0, 0.0])

        r_sq = get_r_sq(y)
        m_paper_rem = get_mass_paper_rem(r_sq)
        M_roll = get_total_mass_roll(r_sq)
        I_roll = get_moment_of_inertia_roll(r_sq, m_paper_rem)
        
        # PE change term N(y) and its derivative M_roll
        N_pe = M_roll * y + lambda_paper * y**2 / 2.0
        # Effective inertia term K(y) and its derivative
        K_ie = M_roll + I_roll / r_sq
        dKdy = get_dK_dy(r_sq)
        
        # Acceleration derived from d(Energy)/dt = 0
        numerator = g * (M_roll * K_ie - N_pe * dKdy)
        denominator = K_ie**2
        
        a = numerator / denominator if denominator else 0.0
        return np.array([v, a])

    # Step 5: Implement 4th Order Runge-Kutta (RK4) integration
    y, v, t = 0.0, 0.0, 0.0
    dt = 0.001  # Time step in seconds
    S = np.array([y, v])
    y_prev, t_prev = y, t

    while y < L:
        y_prev, t_prev = y, t
        k1 = odefunc(t, S)
        k2 = odefunc(t + dt/2, S + dt/2 * k1)
        k3 = odefunc(t + dt/2, S + dt/2 * k2)
        k4 = odefunc(t + dt, S + dt * k3)
        S += (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
        y, v = S[0], S[1]
        t += dt

    # Step 6: Interpolate to find the exact time when y=L for better accuracy
    if y > y_prev:
        t_final = t_prev + dt * (L - y_prev) / (y - y_prev)
    else:
        t_final = t

    # Step 7: Print the final answer and relevant parameters
    print("Based on the provided parameters:")
    print("Diameter of inner cylinder: {} cm".format(d_cylinder * 100))
    print("Toilet paper thickness: {} mm".format(d_paper * 1000))
    print("Number of wraps: {}".format(n_wraps))
    print("Total mass of paper: {} g".format(m_paper_total * 1000))
    print("Mass of cardboard cylinder: {} g".format(m_cylinder * 1000))
    print("\nDerived properties:")
    print("Total length of paper to unroll: {:.2f} m".format(L))
    print("\nResult:")
    print("The time it takes for the toilet paper to reach the end of its roll is: {:.2f} seconds.".format(t_final))
    return t_final

if __name__ == '__main__':
    final_time = solve_toilet_paper_fall()
    # Return answer in specified format
    # <<<1.56>>>
    

solve_toilet_paper_fall()
print(f'<<<{solve_toilet_paper_fall():.2f}>>>')