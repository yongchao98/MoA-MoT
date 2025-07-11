import numpy as np
from scipy import integrate

def calculate_cross_section():
    """
    Calculates the total cross-section for fermion-fermion scattering
    via pseudoscalar exchange for a given set of parameters.
    """
    # --- Parameters ---
    # These are example values. You can change them to fit your specific problem.
    # We use units where c=1 and h_bar=1 (GeV for energy/mass).
    g = 1.0       # Dimensionless coupling constant
    E = 10.0      # Energy of each fermion in the CM frame (in GeV)
    m = 0.5       # Mass of the fermion psi (in GeV)
    M = 1.0       # Mass of the scalar boson phi (in GeV)

    print(f"Calculating total cross-section σ for the following parameters:")
    print(f"Coupling constant g = {g}")
    print(f"Fermion mass m = {m} GeV")
    print(f"Scalar mass M = {M} GeV")
    print(f"CM-frame energy per fermion E = {E} GeV\n")

    # Check for kinematic possibility
    if E < m:
        print("Error: The energy E of the fermion must be greater than its mass m.")
        return

    # --- Integrand function for <|M|^2> ---
    # This function calculates the spin-averaged matrix element squared
    # as a function of x = cos(theta).
    def integrand_M2(x, E, m, M, g):
        # Mandelstam variables in the CM frame
        s = 4 * E**2
        # Squared momentum in the CM frame
        p_squared = E**2 - m**2
        t = -2 * p_squared * (1 - x)
        u = -2 * p_squared * (1 + x)

        # Denominators
        den_t = t - M**2
        den_u = u - M**2

        # To avoid division by zero if t=M^2 or u=M^2 (a resonance case)
        if abs(den_t) < 1e-9 or abs(den_u) < 1e-9:
            # This is a pole, the integral is ill-defined without a width for the mediator.
            # We return a large number or handle it as an error. For simplicity, we can skip this point.
            return np.inf

        # Terms of the squared matrix element
        term_t = t**2 / den_t**2
        term_u = u**2 / den_u**2
        term_interference = - (t * (s + t - 2*m**2)) / (den_t * den_u)

        M2_avg = g**4 * (term_t + term_u + term_interference)
        return M2_avg

    # --- Numerical Integration ---
    # We integrate <|M|^2> with respect to d(cos(theta)) from -1 to 1.
    integral_val, integral_err = integrate.quad(
        integrand_M2,
        -1,  # cos(theta) from -1 (theta=pi)
        1,   # to 1 (theta=0)
        args=(E, m, M, g)
    )

    # --- Total Cross-Section Calculation ---
    # Formula: sigma = (1 / (128 * pi * E^2)) * Integral(<|M|^2> d(cos(theta)))
    prefactor = 1 / (128 * np.pi * E**2)
    sigma = prefactor * integral_val

    # As requested, printing the numbers in the final equation
    print("--- Calculation ---")
    print(f"Integral of <|M|^2> d(cos(theta)) from -1 to 1 = {integral_val:.6e}")
    print(f"Prefactor = 1 / (128 * pi * E^2) = 1 / (128 * {np.pi:.4f} * {E**2:.2f}) = {prefactor:.6e}")
    print(f"σ = Prefactor * Integral")
    print(f"σ = {prefactor:.6e} * {integral_val:.6e} = {sigma:.6e} GeV^-2")
    
    # You can convert this to other units, e.g., barns.
    # 1 GeV^-2 = 0.389379 mb (millibarns)
    sigma_mb = sigma * 0.389379e-3 * 1e9 # in pb (pico barns)
    #print(f"σ = {sigma_mb:.6f} pb")


    # Final answer in the required format
    print(f"\nFinal Answer (in GeV^-2):")
    print(f"<<<{sigma}>>>")


if __name__ == '__main__':
    calculate_cross_section()