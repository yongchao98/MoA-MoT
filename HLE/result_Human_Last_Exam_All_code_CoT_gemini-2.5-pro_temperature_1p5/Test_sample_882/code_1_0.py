import math
from scipy.optimize import brentq

def solve_for_p_and_calculate_result():
    """
    This function solves the game theory problem to find the required value.
    """
    # Problem parameters
    N = 8
    num_players = 3 * N
    num_races = N
    
    # In equilibrium, the payoff for any strategy must be 1/3.
    # The expected payoff for the discrete strategy is E_D(p) = (1 - (1 - p/N)**(3*N)) / (3*p).
    # Setting E_D(p) = 1/3 gives the equation: 1 - (1 - p/N)**(3*N) = p.
    
    # We define a function f(p) = 0 that we need to solve for p.
    # f(p) = 1 - (1 - p/N)**(3*N) - p
    def equilibrium_equation(p, N_val):
        """The equation to solve: 1 - (1 - p/N)**(3N) - p = 0"""
        if p < 0 or p > 1:
            return float('nan') # p must be a probability
        # The equation holds for p=0, but we seek the non-trivial equilibrium.
        if p == 0:
            return 0
        
        term1 = 1.0
        term2 = (1.0 - p / N_val)**(3 * N_val)
        term3 = p
        return term1 - term2 - term3

    # The equation has a trivial root at p=0. We look for the non-trivial root in (0, 1].
    # Using a robust numerical solver, brentq, to find the root.
    # We search in the interval [1e-9, 1.0] to avoid the p=0 root.
    try:
        p_solution = brentq(equilibrium_equation, 1e-9, 1.0, args=(N,))
    except (ValueError, RuntimeError):
        print("Could not find a unique solution for p in the interval (0, 1].")
        return

    # Now calculate the final result
    one_minus_p = 1 - p_solution
    final_value = 10000 * one_minus_p
    result = math.floor(final_value)

    # Output the components of the final equation and the solution
    print("For N = 8, the equilibrium equation is:")
    print(f"1 - (1 - p/{N})**{3*N} = p")
    print("\nNumerically solving this equation gives:")
    print(f"p ≈ {p_solution:.6f}")
    print("\nCalculating the final requested value:")
    print(f"⌊10000 * (1 - p)⌋ = ⌊10000 * (1 - {p_solution:.6f})⌋")
    print(f"                   = ⌊10000 * {one_minus_p:.6f}⌋")
    print(f"                   = ⌊{final_value:.6f}⌋")
    print(f"                   = {result}")

# Run the solver and print the result
solve_for_p_and_calculate_result()