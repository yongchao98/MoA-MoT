import math
from scipy.special import comb
from scipy.optimize import brentq

def solve_game_theory_problem():
    """
    Solves for the mixed strategy equilibrium probability p and computes the final result.
    """
    N = 8
    M = 3 * N - 1

    # Define the utility function for the discrete strategy S_D
    def u_d(p, n_val=N):
        """Calculates the expected utility for a player using the discrete strategy."""
        if p == 0:
            return 1.0
        q = p / n_val
        m_val = 3 * n_val - 1
        return (1 - (1 - q)**(m_val + 1)) / ((m_val + 1) * q)

    # Define the utility function for the continuous strategy S_C
    def u_c(p, n_val=N):
        """Calculates the expected utility for a player using the continuous strategy."""
        m_val = 3 * n_val - 1
        
        # The utility is an expectation over the number of discrete opponents 'l'.
        # E[ N * ((N-1)/N)^l / (M-l+1) ] where l ~ Bin(M, p)
        expected_value = 0
        x = (n_val - 1) / n_val
        
        for l in range(m_val + 1):
            # Binomial probability P(L=l)
            try:
                prob_l = comb(m_val, l, exact=True) * (p**l) * ((1-p)**(m_val - l))
            except ValueError: # Handle potential floating point errors for p=0 or p=1
                if (p == 1.0 and l == m_val) or (p == 0.0 and l == 0):
                    prob_l = 1.0
                else:
                    prob_l = 0.0

            # Payoff given l discrete opponents
            payoff_l = n_val * (x**l) / (m_val - l + 1)
            
            expected_value += prob_l * payoff_l
        return expected_value

    # Define the difference function whose root we want to find
    def diff_func(p):
        """Difference between the two utility functions."""
        return u_d(p) - u_c(p)

    # Solve for p in the interval (0, 1) using Brent's method.
    # We confirmed that diff_func(0) > 0 and diff_func(1) < 0, so a root exists.
    p_star = brentq(diff_func, 1e-9, 1-1e-9)

    # Calculate the final result
    value = 10000 * (1 - p_star)
    result = math.floor(value)

    print(f"For N = {N}, we solve the equation U_D(p) = U_C(p).")
    print(f"The equation is: (1 - (1-p/{N})^{3*N})/(3*p) = Sum_{{l=0..{M}}}[C({M},l)*p^l*(1-p)^{{{M}}-l} * ({N}*((({N}-1)/{N}))^l)/({M}-l+1)]")
    
    ud_val = u_d(p_star)
    uc_val = u_c(p_star)
    print(f"\nFound equilibrium probability p* = {p_star:.6f}")
    print(f"At this equilibrium, the payoffs are equal:")
    print(f"  U_D(p*) = {ud_val:.6f}")
    print(f"  U_C(p*) = {uc_val:.6f}")

    print(f"\nNow we calculate the final value:")
    print(f"1 - p* = {1 - p_star:.6f}")
    print(f"10000 * (1 - p*) = {value:.6f}")
    print(f"The floor of this value is {result}.")

solve_game_theory_problem()