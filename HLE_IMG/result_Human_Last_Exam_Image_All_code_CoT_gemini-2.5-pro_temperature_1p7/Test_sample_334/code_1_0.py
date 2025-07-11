import numpy as np
from scipy.optimize import fsolve

def solve_soc_bec_problem():
    """
    This script solves the SOC-BEC problem by:
    1. Defining the parameters for the missing plot based on logical deduction.
    2. Defining the physical quantities E'(k) and E''(k).
    3. Setting up the equation (m1 + m2)/2 = 0, which simplifies to k*E'' + E' = 0.
    4. Numerically solving for k_0*, the smallest positive root of this equation.
    5. Calculating and printing the final result n_0 * k_R* / k_0*.
    """

    # Step 1: Define parameters for the missing set.
    # From the reasoning explained in the plan:
    # Base parameters (n_0=4) are (delta, Omega, k_R) = (4, 6, 2).
    # Missing set parameters are (delta/2, Omega/2, kR/2) variations.
    # The missing set is (delta, Omega/2, k_R).
    delta_star = 4.0
    Omega_star = 3.0
    kR_star = 2.0
    
    # The index of the base plot
    n_0 = 4

    # Step 2: Define functions for the first and second derivatives of the energy E_-(k)
    def S(k, delta, Omega, kR):
        """ The square root term in the energy dispersion relation. """
        return np.sqrt((2 * k * kR - delta / 2)**2 + (Omega / 2)**2)

    def E_prime(k, delta, Omega, kR):
        """ The group velocity v(k) = dE_-(k)/dk. """
        Sk = S(k, delta, Omega, kR)
        # Avoid division by zero if Sk is ever zero (it won't be for Omega > 0)
        if Sk == 0:
            return 2 * k
        return 2 * k - (2 * kR * (2 * k * kR - delta / 2)) / Sk

    def E_double_prime(k, delta, Omega, kR):
        """ The inverse of the diffusive mass m_2(k) = d^2E_-(k)/dk^2. """
        Sk = S(k, delta, Omega, kR)
        if Sk == 0:
            return 2.0 # Limit case
        return 2 - (kR**2 * Omega**2) / (Sk**3)

    # Step 3: Define the equation to solve for k_0*
    # The condition is (m1 + m2) / 2 = 0 => m1 = -m2
    # k/E'(k) = -1/E''(k) => k * E''(k) + E'(k) = 0
    def equation_to_solve(k):
        # We need to solve for k > 0.
        if k <= 0:
            return np.inf # We are interested in positive k
        val = k * E_double_prime(k, delta_star, Omega_star, kR_star) + E_prime(k, delta_star, Omega_star, kR_star)
        return val

    # Step 4: Find the smallest positive root k_0*
    # From manual checks, the root lies between 1.0 and 1.25.
    # We use a starting guess of 1.2 for the numerical solver.
    initial_guess = 1.2
    k0_star_solution = fsolve(equation_to_solve, initial_guess)
    k0_star = k0_star_solution[0]

    # Step 5: Calculate and print the final result
    result = n_0 * kR_star / k0_star

    print("--- Determined and Calculated Values ---")
    print(f"Base plot index, n_0 = {n_0}")
    print(f"Raman wavevector for the missing set, k_R* = {kR_star}")
    print(f"Smallest positive k where (m1+m2)/2 = 0, k_0* ≈ {k0_star:.4f}")
    print("\n--- Final Calculation ---")
    print(f"{n_0} * {kR_star} / {k0_star:.4f} = {result:.4f}")
    print(f"\nThe value of n_0 * k_R* / k_0* is approximately {result:.2f}")


solve_soc_bec_problem()