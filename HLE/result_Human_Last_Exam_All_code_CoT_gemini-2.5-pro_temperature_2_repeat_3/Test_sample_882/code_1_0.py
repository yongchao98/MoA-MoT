import math
from scipy.optimize import brentq

def solve_for_p():
    """
    Solves the game theory problem to find the value of p and the final requested quantity.
    """
    N = 8
    num_players = 3 * N

    # The equation for p is derived from the equilibrium condition E_D(p) = 1/3,
    # which simplifies to p = 1 - (1 - p/N)^(3*N).
    # We define a function f(p) = 0 to find the root.
    def f(p):
        return p - 1 + (1 - p / N)**(num_players)

    # We solve for p in the interval (0, 1). We avoid the trivial p=0 root by
    # searching in a smaller interval like (0.1, 1.0).
    try {
        p_solution = brentq(f, 0.1, 1.0)
    } except ValueError:
        print("Could not find a unique root in the interval.")
        return

    # Calculate the final value as requested by the problem.
    result = 10000 * (1 - p_solution)
    final_answer = math.floor(result)
    
    print("For N = 8, the equilibrium equation to solve is:")
    # Output each number in the final equation as requested
    print(f"p = 1 - (1 - p/{N})**{num_players}")
    print("\nSolving this equation numerically gives the probability 'p' that a player uses the discrete strategy.")
    # Report p to six significant digits
    print(f"p = {p_solution:.6f}")
    
    print("\nThe final required quantity is floor(10000 * (1-p)).")
    print(f"1 - p = {1 - p_solution:.6f}")
    print(f"10000 * (1 - p) = {result:.4f}")
    print(f"The floor of the result is: {final_answer}")

solve_for_p()