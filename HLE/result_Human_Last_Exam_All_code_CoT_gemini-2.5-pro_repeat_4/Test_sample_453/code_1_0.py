import numpy as np
from scipy.special import genlaguerre
from scipy.integrate import quad
from scipy.constants import m_e, c, alpha, Hartree, physical_constants
from math import factorial

def solve():
    """
    Calculates the second-order relativistic energy shift for an electron
    in a hydrogen atom in the n=3, l=2 state.
    """
    # State quantum numbers
    n = 3
    l = 2

    # We use atomic units for the core calculation (hbar=1, m_e=1, e=1, a_0=1).
    # In these units, the speed of light c is 1/alpha.
    
    # Hydrogen radial wave function in atomic units (Z=1)
    def R_nl(n_q, l_q, r):
        if r <= 0:
            return 0
        # Normalization constant
        try:
            norm = np.sqrt((2.0/n_q)**3 * factorial(n_q - l_q - 1) / (2.0 * n_q * factorial(n_q + l_q)))
        except ValueError:
            # Factorial of negative number
            return 0
        
        # Laguerre polynomial L_{n-l-1}^{2l+1}
        laguerre = genlaguerre(n_q - l_q - 1, 2*l_q + 1)
        # Argument for Laguerre and exponential
        rho = 2.0 * r / n_q
        
        return norm * np.exp(-rho / 2.0) * rho**l_q * laguerre(rho)

    # Unperturbed energy levels in atomic units (Hartrees)
    def E_n(n_q):
        return -1.0 / (2.0 * n_q**2)

    # Integrands for matrix elements <k|V|n> and <k|V^2|n> where V=-1/r
    def integrand_V(r, k, n_q, l_q):
        # Integrand for <k|V|n> = integral[ R_kl(r) * (-1/r) * R_nl(r) * r^2 dr ]
        return R_nl(k, l_q, r) * (-1.0/r) * R_nl(n_q, l_q, r) * r**2

    def integrand_V2(r, k, n_q, l_q):
        # Integrand for <k|V^2|n> = integral[ R_kl(r) * (1/r^2) * R_nl(r) * r^2 dr ]
        return R_nl(k, l_q, r) * (1.0/r**2) * R_nl(n_q, l_q, r) * r**2

    # Calculate matrix elements using numerical integration
    memo_V = {}
    memo_V2 = {}
    def get_V_kn(k, n_q, l_q):
        if (k, n_q, l_q) in memo_V:
            return memo_V[(k, n_q, l_q)]
        result, _ = quad(integrand_V, 0, np.inf, args=(k, n_q, l_q))
        memo_V[(k, n_q, l_q)] = result
        return result

    def get_V2_kn(k, n_q, l_q):
        if (k, n_q, l_q) in memo_V2:
            return memo_V2[(k, n_q, l_q)]
        result, _ = quad(integrand_V2, 0, np.inf, args=(k, n_q, l_q))
        memo_V2[(k, n_q, l_q)] = result
        return result

    # Summation for the second-order energy shift
    delta_E2_hartrees = 0.0
    k_max = 50  # Sum over a range of intermediate states for convergence
    E_n_val = E_n(n)

    for k in range(1, k_max + 1):
        if k == n:
            continue

        E_k_val = E_n(k)
        
        # Calculate matrix elements of V and V^2
        V_kn = get_V_kn(k, n, l)
        V2_kn = get_V2_kn(k, n, l)

        # Calculate the matrix element of the perturbation H' in atomic units.
        # H' = -p^4/(8m^3c^2) = -alpha^2 * p^4 / 8 (since m=1, c=1/alpha)
        # p^2 = 2m(H0-V) = 2(H0-V). M_kn = <k|H'|n> = - (alpha^2/2) * <k|(H0-V)^2|n>
        # For k!=n, this simplifies to:
        # M_kn = - (alpha^2 / 2) * [ -(E_k + E_n)*<k|V|n> + <k|V^2|n> ]
        M_kn = - (alpha**2 / 2.0) * (-(E_k_val + E_n_val) * V_kn + V2_kn)

        # Add term to the sum
        delta_E2_hartrees += abs(M_kn)**2 / (E_n_val - E_k_val)

    # The result delta_E2_hartrees is in atomic units (Hartrees)
    # Convert to Joules to find the dimensionless coefficient C
    delta_E2_joules = delta_E2_hartrees * Hartree

    # Express the result in the form C * m_e * c^2 * alpha^6
    factor = m_e * c**2 * alpha**6
    C = delta_E2_joules / factor

    # Output the final result
    print("The second-order energy shift due to the relativistic kinetic energy correction is given by the formula:")
    print("ΔE = C * m_e * c^2 * α^6")
    print("\nFor the state n=3, l=2, the calculated dimensionless coefficient C is:")
    print(f"C = {C}")
    
    # Print the final equation with the calculated number for C
    print("\nTherefore, the final equation for the energy shift is:")
    print(f"ΔE = {C} * m_e * c^2 * α^6")

solve()
<<<ΔE = -0.0003328591823908846 * m_e * c^2 * α^6>>>