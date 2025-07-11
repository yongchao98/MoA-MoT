import numpy as np
from scipy.optimize import brentq

def solve_game():
    """
    Solves for the mixed strategy equilibrium probability p and computes the final result.
    """
    N = 8
    
    def h(p):
        """
        Represents the difference in expected payoffs, H(p) = E_D(p) - E_S(p).
        The root of this function is the equilibrium probability p.
        """
        # Handle boundary cases to avoid division by zero
        if p < 1e-9:
            # lim p->0 H(p) = 1 - 1/3 = 2/3
            return 2/3
        if p > 1 - 1e-9:
            # lim p->1 H(p) = U_D(1) - U_S(1)
            ud1 = (1/3) * (1 - (1 - 1/N)**(3*N))
            us1 = N * (1 - 1/N)**(3*N-1)
            return ud1 - us1

        # Expected payoff for Discrete strategy
        e_d = (1 - (1 - p/N)**(3*N)) / (3 * p)
        
        # Expected payoff for Split strategy
        term_s1 = (1 - p/N)**(3*N)
        term_s2 = p**(3*N) * (1 - 1/N)**(3*N)
        e_s = (term_s1 - term_s2) / (3 * (1 - p))
        
        return e_d - e_s

    # Find the root of h(p) in the interval (0, 1) to find the equilibrium p
    # We know a root exists because h(p) is continuous, positive near 0, and negative near 1.
    p_star = brentq(h, a=0, b=1)

    # The problem asks for floor(10000 * (1-p))
    one_minus_p = 1 - p_star
    multiplied_value = 10000 * one_minus_p
    final_answer = np.floor(multiplied_value)

    print(f"For N = {N}, the equation for the equilibrium probability p is E_D(p) = E_S(p), where:")
    print(f"E_D(p) = (1 - (1 - p/{N})^{3*N}) / (3*p)")
    print(f"E_S(p) = ((1 - p/{N})^{3*N} - p^{3*N} * (1 - 1/{N})^{3*N}) / (3*(1-p))")
    print(f"\nSolving this numerically gives p = {p_star:.10f}")
    
    print("\nNow we calculate the final answer as floor(10000 * (1-p)):")
    print(f"1 - p = 1 - {p_star:.10f} = {one_minus_p:.10f}")
    print(f"10000 * (1 - p) = 10000 * {one_minus_p:.10f} = {multiplied_value:.10f}")
    print(f"floor({multiplied_value:.10f}) = {final_answer}")

solve_game()