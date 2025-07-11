import math
from scipy.optimize import root_scalar

def solve_game():
    """
    Solves the game theory problem to find the equilibrium probability p
    and computes the final requested value.
    """
    N = 8
    n_players = 3 * N

    def f(k, p):
        """
        Helper function derived from the game's payoff formula.
        """
        # The limit of f(k, p) as k->0 or p->0 is 1.
        if k == 0 or p == 0:
            return 1.0
        
        base = 1 - p * k / N
        numerator = 1 - base**n_players
        denominator = 3 * p * k
        return numerator / denominator

    def equation_to_solve(p):
        """
        This function represents the equilibrium condition pi(S_N) = pi(D).
        We are looking for the root of this function for p in (0, 1).
        The equation is: sum_{j=1 to N} (-1)^{j-1} * C(N,j) * f(N-j, p) - f(1, p) = 0
        For N=8, this expands to:
        8*f(7,p) - 28*f(6,p) + 56*f(5,p) - 70*f(4,p) + 56*f(3,p) - 28*f(2,p) + 8*f(1,p) - f(0,p) - f(1,p) = 0
        which simplifies to:
        8*f(7,p) - 28*f(6,p) + 56*f(5,p) - 70*f(4,p) + 56*f(3,p) - 28*f(2,p) + 7*f(1,p) - 1 = 0
        """
        if p <= 0 or p >= 1:
            return float('nan')
            
        val = (
            math.comb(N, 1) * f(N - 1, p)
            - math.comb(N, 2) * f(N - 2, p)
            + math.comb(N, 3) * f(N - 3, p)
            - math.comb(N, 4) * f(N - 4, p)
            + math.comb(N, 5) * f(N - 5, p)
            - math.comb(N, 6) * f(N - 6, p)
            + math.comb(N, 7) * f(N - 7, p)
            - math.comb(N, 8) * f(N - 8, p)
            - f(1, p)
        )
        return val

    # We solve the equation for p in the interval (0, 1).
    # Analysis shows that equation_to_solve(epsilon) is negative for small epsilon > 0,
    # and equation_to_solve(1) is positive. So a root exists in (0, 1).
    # We use a robust numerical solver to find it.
    try:
        sol = root_scalar(equation_to_solve, bracket=[0.01, 0.999], method='brentq')
        p = sol.root
    except (ValueError, RuntimeError) as e:
        print(f"Solver failed: {e}")
        return

    # Calculate the final result as requested by the problem.
    result = math.floor(10000 * (1 - p))
    
    print(result)

solve_game()