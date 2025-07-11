import numpy as np
from scipy import constants, integrate

def solve_photon_gas_equilibrium():
    """
    Calculates and prints the equilibrium mean energy and entropy of a photon gas
    in a given volume and at a given temperature.
    
    The calculation is based on the principles of statistical mechanics for a
    Bose gas (photons) and uses numerical integration over the density of states.
    """

    # --- Define System Parameters ---
    T = 500.0  # Temperature in Kelvin (e.g., surface of a hot plate)
    V = 1.0    # Volume in cubic meters

    print(f"--- System Parameters ---")
    print(f"Temperature (T): {T} K")
    print(f"Volume (V):      {V} m^3\n")

    # --- Physical Constants ---
    k_B = constants.k      # Boltzmann constant (J/K)
    hbar = constants.hbar  # Reduced Planck constant (J*s)
    c = constants.c        # Speed of light (m/s)

    # Lagrange multiplier beta
    beta = 1.0 / (k_B * T)

    # --- Density of States Function ---
    # g(epsilon) = (V / (pi^2 * hbar^3 * c^3)) * epsilon^2
    # We define the constant factor first for efficiency.
    g_factor = V / (np.pi**2 * hbar**3 * c**3)

    # --- Integrands for Energy and Entropy Calculation ---
    
    # Integrand for mean energy E
    # integrand_E = epsilon * g(epsilon) / (exp(beta*epsilon) - 1)
    def integrand_E(epsilon):
        # Use np.expm1 for better numerical stability for small arguments
        return (g_factor * epsilon**3) / np.expm1(beta * epsilon)

    # Integrand for the log(Z) term in the entropy calculation
    # log(Z) = -Integral[g(epsilon) * log(1 - exp(-beta*epsilon))] d_epsilon
    # So we integrate g(epsilon) * log(1 - exp(-beta*epsilon))
    def integrand_logZ_term(epsilon):
        # Use np.log1p for better numerical stability
        return g_factor * epsilon**2 * np.log1p(-np.exp(-beta * epsilon))

    # --- Numerical Integration ---
    # We integrate from 0 to infinity
    E, E_err = integrate.quad(integrand_E, 0, np.inf)
    
    # Note: log(Z) is negative. The integral computes the part without the leading minus sign.
    logZ_term, logZ_err = integrate.quad(integrand_logZ_term, 0, np.inf)

    # The total partition function log is log(Z) = -logZ_term
    log_Z = -logZ_term
    
    # --- Calculate Final Equilibrium Values ---
    # S = E/T + k_B * log(Z)
    S = (E / T) + (k_B * log_Z)

    # --- Print Results ---
    print("--- Equilibrium Values ---")
    print(f"Mean Energy (E): {E:.4e} Joules")
    print(f"Entropy (S):     {S:.4e} J/K\n")
    
    # The final equation for Entropy is S = E/T + k_B * log(Z)
    print("--- Entropy Equation Breakdown (S = E/T + k_B*log(Z)) ---")
    term1_E_over_T = E / T
    term2_kBlogZ = k_B * log_Z
    
    print(f"Term 1 (E/T):        {term1_E_over_T:.4e} J/K")
    print(f"Term 2 (k_B*log(Z)):  {term2_kBlogZ:.4e} J/K")
    print(f"Calculated Sum (S):  {term1_E_over_T + term2_kBlogZ:.4e} J/K")

    # The Stefan-Boltzmann Law states E/V = a*T^4. Let's verify our constant 'a'.
    a_radiation_constant = (np.pi**2 * k_B**4) / (15 * hbar**3 * c**3)
    E_analytical = a_radiation_constant * V * T**4
    print(f"\n--- Analytical Check (Stefan-Boltzmann Law) ---")
    print(f"Numerical E:   {E:.4e} J")
    print(f"Analytical E:  {E_analytical:.4e} J")

    return E, S

if __name__ == '__main__':
    E_final, S_final = solve_photon_gas_equilibrium()
    # Output the final answer in the required format
    print(f"\n<<<{E_final:.4e}, {S_final:.4e}>>>")
