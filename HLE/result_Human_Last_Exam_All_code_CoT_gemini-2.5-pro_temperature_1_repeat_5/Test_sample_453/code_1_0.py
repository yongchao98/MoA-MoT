import numpy as np
from scipy.special import genlaguerre
from scipy.integrate import quad
from math import factorial

def solve_energy_shift():
    """
    Calculates the second-order energy shift for an electron in a hydrogen atom
    due to the relativistic kinetic energy correction.
    """
    # Use atomic units where hbar = m_e = e = 1/(4*pi*epsilon_0) = 1.
    # In these units, the speed of light c is the reciprocal of the fine-structure constant.
    alpha = 1 / 137.035999084
    c = 1 / alpha

    # Define the state of interest.
    n_initial = 3
    l_initial = 2

    # The unperturbed energy levels of the hydrogen atom in atomic units.
    def E0(n):
        """Calculates the unperturbed energy E_n^{(0)} = -1/(2n^2) in Hartrees."""
        return -1.0 / (2 * n**2)

    # The radial wavefunction R_nl(r) for the hydrogen atom in atomic units.
    def R_nl(n, l, r):
        """
        Calculates the value of the radial wavefunction R_nl at radius r.
        """
        if n <= l:
            return 0
        norm = np.sqrt((2.0 / n)**3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l)))
        rho = 2.0 * r / n
        laguerre_poly = genlaguerre(n - l - 1, 2 * l + 1)
        return norm * np.exp(-rho / 2.0) * (rho**l) * laguerre_poly(rho)

    # The perturbation H' can be expressed as -1/(2*m*c^2) * (H_0 - V)^2.
    # In atomic units (m=1), the matrix element <n'l'|H'|nl> for n' != n is:
    # <n'l'|H'|nl> = (1/(2c^2)) * [ (E_n' + E_n)<n'l'|V|nl> - <n'l'|V^2|nl> ]
    # Since V(r) is a scalar operator, matrix elements are non-zero only for l'=l.
    
    memo_h_prime = {}
    def calculate_h_prime_matrix_element(n_prime, l_prime, n, l):
        """
        Calculates the off-diagonal matrix element <n',l'|H'|n,l>.
        """
        if (n_prime, l_prime, n, l) in memo_h_prime:
            return memo_h_prime[(n_prime, l_prime, n, l)]

        if l_prime != l or n_prime == n:
            return 0

        # V = -1/r, V^2 = 1/r^2 in atomic units.
        V_integrand = lambda r: R_nl(n_prime, l, r) * (-1.0 / r) * R_nl(n, l, r) * r**2
        V2_integrand = lambda r: R_nl(n_prime, l, r) * (1.0 / r**2) * R_nl(n, l, r) * r**2

        V_me, _ = quad(V_integrand, 0, np.inf, limit=200)
        V2_me, _ = quad(V2_integrand, 0, np.inf, limit=200)
        
        E_n_prime = E0(n_prime)
        E_n = E0(n)

        h_prime_me = (1.0 / (2.0 * c**2)) * ((E_n_prime + E_n) * V_me - V2_me)
        memo_h_prime[(n_prime, l_prime, n, l)] = h_prime_me
        return h_prime_me

    # Calculate the second-order energy shift E^(2).
    # E^(2) = sum_{n' != n, l'=l} |<n'l|H'|nl>|^2 / (E_n - E_n')
    # We approximate by summing over a finite number of discrete states.
    E_initial = E0(n_initial)
    E2_shift = 0.0
    N_max_cutoff = 40

    print(f"Calculating the second-order energy shift for the |n={n_initial}, l={l_initial}> state.")
    print(f"Approximating by summing over intermediate states |n', l={l_initial}> up to n'={N_max_cutoff}.")

    for n_prime in range(l_initial + 1, N_max_cutoff + 1):
        if n_prime == n_initial:
            continue
        
        E_intermediate = E0(n_prime)
        me = calculate_h_prime_matrix_element(n_prime, l_initial, n_initial, l_initial)
        term = (me**2) / (E_initial - E_intermediate)
        E2_shift += term

    # The result is in atomic units (Hartrees). The expected form is C * m * c^2 * alpha^6.
    # In atomic units, this energy is C * (1/alpha^2) * alpha^6 = C * alpha^4.
    # So, the numerical coefficient C = E2_shift / alpha^4.
    numerical_coefficient = E2_shift / alpha**4

    print("\nCalculation complete.")
    print(f"The calculated shift in atomic units (Hartrees) is: {E2_shift:.6e}")
    print("The energy shift in terms of fundamental constants is expressed as E^(2) = C * m_e * c^2 * alpha^6.")
    print(f"The calculated numerical coefficient is C = {numerical_coefficient:.5f}")
    
    print("\nThe final expression for the energy shift is:")
    print(f"E^(2) = {numerical_coefficient:.5f} * m_e * c^2 * \u03b1^6")
    
    # Returning the coefficient for the final answer tag.
    return numerical_coefficient

# Run the calculation and store the coefficient for the final answer.
final_coefficient = solve_energy_shift()