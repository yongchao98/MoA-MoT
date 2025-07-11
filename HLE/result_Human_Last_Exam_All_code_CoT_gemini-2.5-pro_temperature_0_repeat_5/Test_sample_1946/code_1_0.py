import numpy as np
from scipy.integrate import quad

def solve_toilet_paper_problem():
    """
    Calculates the time it takes for a falling toilet paper roll to completely unroll.
    
    The solution is based on the principle of energy conservation and uses numerical
    integration to find the total time.
    """
    
    # 1. Define physical constants from the problem statement.
    g = 9.81  # Acceleration due to gravity (m/s^2)
    d_inner = 0.04  # Inner cylinder diameter (m)
    r_inner = d_inner / 2.0  # Inner cylinder radius (m)
    t_paper = 0.0005  # Paper thickness (m)
    N_wraps = 100  # Number of wraps
    m_paper = 0.2  # Total mass of paper (kg)
    m_cylinder = 0.02  # Mass of cardboard cylinder (kg)

    # 2. Calculate derived constants.
    # Initial outer radius of the full roll.
    R0 = r_inner + N_wraps * t_paper
    
    # Total length of the paper, calculated from the volume.
    # Volume = Area_of_annulus * width = L_total * t_paper * width
    # L_total = Area_of_annulus / t_paper
    L_total = np.pi * (R0**2 - r_inner**2) / t_paper
    
    # Linear mass density of the paper.
    rho_L = m_paper / L_total
    
    # Moment of inertia of the cardboard cylinder, assuming it's a thin shell (I = mr^2).
    I_cylinder = m_cylinder * r_inner**2

    # Print the parameters used in the calculation.
    print("--- Parameters Used in Calculation ---")
    print(f"Inner cylinder diameter: {d_inner * 100:.1f} cm")
    print(f"Paper thickness: {t_paper * 1000:.1f} mm")
    print(f"Number of wraps: {N_wraps}")
    print(f"Mass of paper: {m_paper * 1000:.0f} grams")
    print(f"Mass of cylinder: {m_cylinder * 1000:.0f} grams")
    print(f"Initial outer radius: {R0 * 100:.1f} cm")
    print(f"Total paper length: {L_total:.2f} m")
    print("------------------------------------")

    # 3. Define functions for the roll's properties as a function of unrolled length y.
    def get_mass_on_roll(y):
        """Mass of the paper remaining on the roll."""
        if y >= L_total:
            return 0
        return m_paper * (1 - y / L_total)

    def get_total_mass(y):
        """Total mass of the falling object (cylinder + paper on roll)."""
        return m_cylinder + get_mass_on_roll(y)

    def get_radius_squared(y):
        """Square of the current outer radius of the roll."""
        if y >= L_total:
            return r_inner**2
        # This is derived from the ratio of remaining paper area to total paper area.
        return r_inner**2 + (R0**2 - r_inner**2) * (1 - y / L_total)

    def get_moment_of_inertia(y):
        """Total moment of inertia of the falling object."""
        if y >= L_total:
            return I_cylinder
        m_p_roll = get_mass_on_roll(y)
        R_sq = get_radius_squared(y)
        # Moment of inertia of the paper on the roll (a hollow cylinder).
        I_p_roll = 0.5 * m_p_roll * (R_sq + r_inner**2)
        return I_cylinder + I_p_roll

    # 4. Define the integrand for the time calculation, which is 1/v(y).
    def integrand(y):
        """
        Calculates the inverse of the velocity, 1/v(y).
        v(y) is derived from the energy conservation equation.
        """
        # The case y=0 is a singularity (v=0), but the integral is convergent.
        # The quad function can handle this type of singularity.
        if y <= 0:
            return np.inf

        M_y = get_total_mass(y)
        R_sq_y = get_radius_squared(y)
        I_y = get_moment_of_inertia(y)

        # Numerator of v^2, from the change in potential energy.
        # PE_roll + PE_unrolled_paper
        pe_term = M_y * g * y + 0.5 * rho_L * g * y**2
        
        # Denominator of v^2, from the kinetic energy terms.
        # KE_trans + KE_rot + KE_unrolled_paper
        ke_coeff = 0.5 * M_y + 0.5 * I_y / R_sq_y + (1/6.0) * rho_L * y

        if ke_coeff <= 0 or pe_term < 0:
            return np.inf

        v_squared = pe_term / ke_coeff
        v = np.sqrt(v_squared)
        
        return 1.0 / v

    # 5. Perform the numerical integration to find the total time.
    # We integrate dt = dy/v(y) from y=0 to y=L_total.
    time, error = quad(integrand, 0, L_total)

    print(f"\nIt takes {time:.2f} seconds for the toilet paper to unroll completely.")
    
    return time

# Execute the solver and store the final answer.
final_time = solve_toilet_paper_problem()
print(f"\n<<< {final_time:.2f} >>>")
