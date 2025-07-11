import numpy as np
from scipy.optimize import fsolve

def solve_game_equilibrium():
    """
    Solves for the mixed strategy equilibrium probability p and computes the final value.
    """
    N = 8
    TOTAL_PLAYERS = 3 * N

    # Define the equilibrium equation f(p) = 0 based on the indifference principle.
    # (1-p)(1 - (1-p/N)^3N) - p((1-p/N)^3N - (1-(2-p)/N)^3N) = 0
    def equilibrium_equation(p):
        # We are looking for a solution where p is a probability between 0 and 1.
        if p <= 0 or p >= 1:
            return 1e9 # Return a large value to guide solver away from boundaries

        term1 = (1 - p / N) ** TOTAL_PLAYERS
        term2 = (1 - (2 - p) / N) ** TOTAL_PLAYERS
        
        # The equation derived from E_D(p) = E_S(p)
        equation = (1 - p) * (1 - term1) - p * (term1 - term2)
        return equation

    # As established in the literature, for N=8 a mixed strategy equilibrium exists.
    # We need to find the root of the equilibrium equation.
    # We use a numerical solver, starting with an initial guess.
    # The solution is known to be close to 1. Let's start the search at 0.95.
    initial_guess = 0.95
    p_solution = fsolve(equilibrium_equation, initial_guess)[0]

    # Calculate 1-p
    one_minus_p = 1 - p_solution
    
    # Calculate the final result
    result = np.floor(10000 * one_minus_p)

    print("For N = 8:")
    print(f"The equilibrium probability p is approximately: {p_solution:.6f}")
    
    # The final equation is floor(10000 * (1 - p))
    # We need to show the components of the final calculation
    print(f"The value of 1-p is: {one_minus_p:.6f}")
    print(f"The value of 10000 * (1-p) is: {10000 * one_minus_p:.6f}")
    print(f"The final result is floor(10000 * (1-p)) = {int(result)}")

solve_game_equilibrium()