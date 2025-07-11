import math
from scipy.optimize import root_scalar

def solve_game_theory_problem():
    """
    This function solves for the probability 'p' in the mixed strategy Nash Equilibrium
    and computes the final result as requested by the user.
    """
    N = 8
    n_players = 3 * N
    
    # Payoff for the discrete strategy (D)
    def u_d(p):
        if p == 0:
            # As p -> 0, the limit of the payoff is 1.
            return 1.0
        # For p > 0, use the formula
        return (1 - (1 - p / N)**n_players) / (3 * p)

    # Payoff for the strategy of splitting fuel among m races (Sm)
    # For N=8, the optimal deviation is m=8.
    def u_s8(p):
        m = 8
        if p == 0:
            # As p -> 0, the limit is 1.
            return 1.0
        
        payoff = 0.0
        for j in range(1, m + 1):
            try:
                # Binomial coefficient C(m, j)
                comb = math.comb(m, j)
                # The term (1 - j*p/N)
                base = 1 - j * p / N
                # The full term in the summation
                term = ((-1)**(j - 1)) * comb * (base**(n_players - 1))
                payoff += term
            except ValueError:
                # This can happen if the base becomes negative, which is possible for large j and p
                # In this specific problem context (p near 1), it shouldn't be an issue.
                pass
        return payoff

    # We need to find the root of the function f(p) = U_D(p) - U_S8(p)
    def f(p):
        return u_d(p) - u_s8(p)

    # Based on analysis of the problem, the solution for p is expected to be close to 1.
    # We will use a numerical root-finding algorithm. `root_scalar` is a robust choice.
    # It requires a bracket [a, b] where f(a) and f(b) have opposite signs.
    # Analysis shows f(p) is negative for p near 0 and at p=1, but a solution must exist.
    # This suggests the function becomes positive somewhere in between.
    # Let's search for a sign change to create a bracket for the solver.
    # If a simple bracket like [0.01, 1.0] fails, it indicates a more complex function shape.
    # However, for this problem, the root is known to be in this interval.
    try:
        # We search for the non-trivial root in the interval (0, 1)
        sol = root_scalar(f, bracket=[0.01, 1.0], method='brentq')
        p = sol.root
    except ValueError:
        print("Solver failed to find a root within the bracket. The function may not cross zero as expected.")
        print("This may indicate a subtle issue in the problem model or formulas.")
        # As a fallback, use a method that doesn't require a bracket, like Newton's method.
        # This is less robust and needs a good starting guess.
        try:
            sol = root_scalar(f, x0=0.9, method='newton')
            p = sol.root
        except (RuntimeError, ValueError) as e:
            print(f"Newton's method also failed: {e}")
            return

    # Calculate the final result
    result = math.floor(10000 * (1 - p))
    
    # Output the final equation with the calculated value of p
    ud_val = u_d(p)
    us8_val = u_s8(p)
    print(f"For N=8, the equilibrium condition is U_D(p) = U_S8(p).")
    print(f"We solve the equation numerically:")
    print(f"(1 - (1 - p/{N})**{n_players}) / (3*p) = Sum_{{j=1 to 8}} [(-1)**(j-1) * C(8,j) * (1 - j*p/{N})**{n_players-1}]")
    print(f"\nThe solution is p = {p:.6f}")
    print(f"At this value of p:")
    print(f"U_D(p)  = {ud_val:.6f}")
    print(f"U_S8(p) = {us8_val:.6f}")
    print("\nCalculating the final answer:")
    print(f"floor(10000 * (1 - {p:.6f})) = floor({10000 * (1 - p):.6f}) = {result}")

solve_game_theory_problem()