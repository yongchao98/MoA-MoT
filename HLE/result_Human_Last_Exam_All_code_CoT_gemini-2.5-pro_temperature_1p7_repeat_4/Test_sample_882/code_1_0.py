import numpy as np
from scipy.optimize import brentq

def solve_for_p():
    """
    Solves the game theory problem to find the value of p and the final expression.
    """
    N = 8
    M = 3 * N  # Total number of players is 24

    # The equilibrium equation is 1 - (1 - p/N)**(M) = p.
    # We define a function whose root we want to find.
    # Let x = p/N, so p = N*x. The equation becomes:
    # 1 - (1 - x)**M - N*x = 0
    # 1 - (1 - x)**24 - 8*x = 0
    def equation(x):
        return 1 - (1 - x)**M - N * x

    # We need to find the non-trivial root of this equation in the interval (0, 1/N).
    # x=0 is a trivial solution. The other root lies between 0 and 1/N.
    # We use the brentq method for robust root finding.
    try:
        # Search for the root in the interval [1e-9, 1/N].
        x = brentq(equation, 1e-9, 1/N)
    except (ValueError, RuntimeError):
        print("Root finding failed. Could not solve for x.")
        return

    # The probability p is N*x
    p = N * x
    
    # The problem asks for the value of floor(10000 * (1-p))
    one_minus_p = 1 - p
    value_to_floor = 10000 * one_minus_p
    result = np.floor(value_to_floor)

    # As requested, output the numbers in the final equation
    print(f"For N = {N}, the equilibrium equation is: 1 - (1 - p/{N})^{{{M}}} = p")
    print(f"Solving this numerically gives a value for p.")
    print(f"p = {p:.8f}")
    print(f"1 - p = {one_minus_p:.8f}")
    print(f"10000 * (1 - p) = {value_to_floor:.8f}")
    print(f"The final result floor(10000 * (1 - p)) is:")
    print(int(result))

solve_for_p()