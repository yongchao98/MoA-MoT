import numpy as np
from scipy.optimize import brentq

def solve_blotto_game():
    """
    Solves the game theory problem to find the value of floor(10000 * (1-p)).
    """
    N = 8
    
    # These are the coefficients in the equation derived from the model.
    n = 3 * N - 1
    m = 6 * N - 1

    print(f"The problem reduces to solving the following equation for x, where p = {N}*x:")
    print(f"(1 - x)^({n}) = 1 / (1 + ({m})*x)\n")

    def equation_to_solve(x):
      """The transcendental equation that defines the equilibrium condition."""
      # This function is designed to find the root of f(x) = 0.
      return (1 - x)**n - 1 / (1 + m * x)

    # We need to find a non-zero solution for x in the interval (0, 1/N).
    # A root is guaranteed to exist in this interval because the function value
    # changes sign from positive near x=0 to negative at x=1/N.
    try:
        # brentq is a robust root-finding algorithm for a bracketed root.
        x_solution = brentq(equation_to_solve, 1e-9, 1/N)
        
        # Calculate p from the solution for x.
        p = N * x_solution
        
        # Calculate the final value as requested.
        final_value = np.floor(10000 * (1 - p))
        
        print(f"Numerical solution for x: {x_solution:.6f}")
        print(f"The probability p is: {p:.6f}")
        print(f"The value of 1-p is: {1-p:.6f}")
        print(f"The value of 10000 * (1-p) is: {10000 * (1-p):.4f}")
        print(f"The final answer is the floor of this value.")
        print(f"\nFinal Answer: {int(final_value)}")

    except (ValueError, RuntimeError) as e:
        print(f"An error occurred during the numerical solution: {e}")

solve_blotto_game()