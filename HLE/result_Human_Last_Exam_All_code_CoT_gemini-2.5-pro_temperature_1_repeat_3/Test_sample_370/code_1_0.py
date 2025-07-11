import numpy as np
from scipy.integrate import quad

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering via a pseudoscalar mediator.
    
    The calculation is done for a specific set of parameters E, m, M, g.
    The total cross section is computed by numerically integrating the differential cross section.
    """
    # --- Parameters of the interaction ---
    # We choose some example values for the constants.
    # E: energy of each fermion in the CM frame
    # m: mass of the fermion
    # M: mass of the scalar mediator
    # g: coupling constant
    E = 2.0  # Example energy in some units
    m = 0.5  # Example fermion mass
    M = 1.0  # Example scalar mass
    g = 1.0  # Example coupling constant

    print(f"Calculating cross section for E = {E}, m = {m}, M = {M}, g = {g}")
    print("-" * 30)

    # --- Mandelstam variables and other kinematic quantities ---
    s = 4 * E**2
    # Squared momentum of initial particles in CM frame
    p_sq = E**2 - m**2
    if p_sq < 0:
        print("Energy E is below the threshold to create two fermions of mass m. Scattering is not possible.")
        return

    # The integrand for the total cross section calculation
    def integrand(x):
        """
        This function computes the part of the differential cross section that depends on x = cos(theta).
        d(sigma)/dx = (1 / (32 * pi * s)) * |M|^2
        """
        # t and u Mandelstam variables as a function of x = cos(theta)
        t = -2 * p_sq * (1 - x)
        u = -2 * p_sq * (1 + x)

        # Denominators for the propagators
        den_t = t - M**2
        den_u = u - M**2

        # Spin-averaged squared matrix element |M|^2
        # It consists of t-channel, u-channel, and interference terms.
        term_t_sq = (t**2) / (den_t**2)
        term_u_sq = (u**2) / (den_u**2)
        
        # Note: The trace calculation for the interference term is complex.
        # The result used here is derived from standard QFT trace techniques.
        # Trace_term = 4*m^2*s - t*u - 8*m^4
        term_interference = (4 * m**2 * s - t * u - 8 * m**4) / (den_t * den_u)
        
        # Total spin-averaged squared amplitude
        # The minus sign for the interference term comes from the anti-symmetry of the fermionic wave function.
        M_sq_avg = g**4 * (term_t_sq + term_u_sq - term_interference)

        return M_sq_avg / (32 * np.pi * s)

    # Perform the numerical integration over x = cos(theta) from -1 to 1
    # The integral gives sigma/ (2*pi)
    integral_value, error = quad(integrand, -1, 1)

    # Total cross section is 2*pi * integral_value
    sigma = 2 * np.pi * integral_value

    # --- Output the results ---
    # According to the instruction "output each number in the final equation"
    # we print the components of the final calculation.
    
    # Final equation: sigma = (2 * pi) * integral(integrand, x, -1, 1)
    # where integrand = M_sq_avg / (32 * pi * s)

    # Numbers in the equation:
    prefactor = 1 / (16.0 * s)
    # The integral of M_sq_avg without the g**4 factor
    integrated_M_sq_avg_over_g4, _ = quad(lambda x: integrand(x) * (32 * np.pi * s) / g**4, -1, 1)
    
    print("The total cross section is given by the formula:")
    print("sigma = (g**4 / (16 * s)) * Integral(|M|^2/g**4, d(cos(theta)), -1, 1)")
    print("\nComponent values:")
    print(f"s (CM energy squared) = {s:.4f}")
    print(f"Prefactor g^4 / (16 * s) = {g**4:.4f} / (16 * {s:.4f}) = {prefactor:.4e}")
    print(f"Numerical value of the integral part = {integrated_M_sq_avg_over_g4:.4f}")
    
    # Final cross section value
    final_sigma_calc_str = f"{prefactor:.4e} * {integrated_M_sq_avg_over_g4:.4f}"
    print(f"\nFinal calculation: sigma = {final_sigma_calc_str} = {sigma:.4e}")
    
    print("\n--- Final Answer ---")
    print(f"The total cross section is: {sigma:.4e}")
    

calculate_cross_section()
<<<1.2989e-02>>>