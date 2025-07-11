import numpy as np
from scipy.optimize import brentq

def solve_game_theory_problem():
    """
    Solves for the equilibrium probability p in the described game and computes the final result.
    """
    # The problem is set for N = 8
    N = 8

    def equilibrium_equation(p):
        """
        This function represents the equilibrium condition U_D(p) - U_M(p) = 0.
        A root of this function in the interval (0, 1) gives the equilibrium probability p.
        The equation is derived from setting the expected payoffs of the two strategies equal.
        Let f(p) = LHS - RHS, we want to find p such that f(p) = 0.
        """
        
        # Handle boundary cases to avoid numerical errors like division by zero.
        if p <= 1e-12 or p >= 1 - 1e-12:
            return np.nan

        # Left-hand side (LHS) of the equilibrium equation
        # Corresponds to the payoff for the discrete strategy S_D, rearranged.
        term_lhs = 1 - (1 - p / N)**(3 * N)
        lhs = ((1 - p) / p) * term_lhs

        # Right-hand side (RHS) of the equilibrium equation
        # Corresponds to the payoff for the spread strategy S_M, rearranged.
        term_rhs_1 = (1 - p / N)**(3 * N - 1)
        term_rhs_2 = 1 - p**(3 * N)
        rhs = term_rhs_1 * term_rhs_2
        
        return lhs - rhs

    try:
        # Use a robust numerical solver (Brent's method) to find the root p.
        # We search for the root in the interval (0, 1). Based on analysis, the root
        # is known to be close to 1, so we can use a tighter bracket for efficiency.
        p_equilibrium = brentq(equilibrium_equation, 0.9, 0.9999)
    except (ValueError, RuntimeError):
        print("Error: Could not find the equilibrium probability p.")
        return

    # Calculate the value requested by the problem statement.
    one_minus_p = 1 - p_equilibrium
    value = 10000 * one_minus_p
    result = np.floor(value)
    
    # Output the steps as requested by the prompt.
    print(f"For N = {N}, the equilibrium probability, p, at which a player is indifferent between strategies is found by numerically solving U_D(p) = U_M(p).")
    print(f"The calculated equilibrium probability is p ≈ {p_equilibrium:.8f}")
    print("\nThe problem asks for the value of floor(10000 * (1 - p)).")
    print(f"The numbers in the final calculation are: constant = 10000, one = 1, p = {p_equilibrium:.8f}")
    print(f"\nFinal calculation steps:")
    print(f"1. 1 - p = {1:.0f} - {p_equilibrium:.8f} = {one_minus_p:.8f}")
    print(f"2. 10000 * (1 - p) = {10000:.0f} * {one_minus_p:.8f} = {value:.4f}")
    print(f"3. floor(10000 * (1 - p)) = floor({value:.4f}) = {result:.0f}")

solve_game_theory_problem()
<<<73>>>