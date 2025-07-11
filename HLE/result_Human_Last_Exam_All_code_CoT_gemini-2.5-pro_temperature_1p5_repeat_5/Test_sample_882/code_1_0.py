import math
import scipy.optimize

def solve_for_p_and_calculate_value():
    """
    This script solves the equilibrium equation for the described game to find p,
    and then calculates the final requested value.
    """
    
    # Game parameters
    N = 8
    M = 3 * N

    # The equilibrium equation is: (1 - p/N)^M = 1 - p
    # We define a function f(p) = (1 - p/N)^M - (1 - p) whose root we need to find.
    def f(p):
        return (1 - p / N)**M - (1 - p)

    # We need to find the non-trivial root of f(p) = 0 in the interval (0, 1).
    # We use a numerical solver, scipy.optimize.brentq, which is robust for this task.
    # We search in the interval [0.1, 1.0] as p=0 is a trivial solution.
    try:
        p_solution = scipy.optimize.brentq(f, 0.1, 1.0)
    except ValueError:
        print("Could not find a root in the specified interval.")
        return

    # The final step is to calculate the expression floor(10000 * (1-p)).
    one_minus_p = 1 - p_solution
    final_value = math.floor(10000 * one_minus_p)

    # Output the results as requested.
    print(f"The equation to solve is (1 - p/N)^M = 1 - p")
    print(f"With N = {N} and M = {M}, this becomes (1 - p/{N})^{M} = 1 - p.")
    print(f"The probability 'p' is solved numerically.")
    print(f"p = {p_solution:.6f}")
    print(f"1 - p = {one_minus_p:.6f}")
    print(f"10000 * (1 - p) = {10000 * one_minus_p:.6f}")
    print(f"The final result floor(10000 * (1 - p)) is: {final_value}")

solve_for_p_and_calculate_value()