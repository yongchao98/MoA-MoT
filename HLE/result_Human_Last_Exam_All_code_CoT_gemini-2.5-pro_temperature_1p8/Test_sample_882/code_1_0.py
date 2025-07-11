import math
from scipy.optimize import root_scalar

def solve_for_p_and_calculate():
    """
    Solves for the mixed strategy equilibrium probability p and calculates the final value.
    """
    # Set the number of races
    N = 8
    
    # Define the equation for which we need to find the root.
    # The equation is derived from the equilibrium condition U_D(p) = 1/3,
    # which simplifies to p = 1 - (1 - p/N)^(3N), or f(p) = 0 where:
    # f(p) = p + (1 - p/N)^(3*N) - 1
    def equation_to_solve(p):
        return p + (1 - p / N)**(3 * N) - 1

    # Solve the equation numerically for p. We search for a root in the interval (0, 1).
    # Initial manual calculations suggest the root is near 0.95.
    try:
        sol = root_scalar(equation_to_solve, bracket=[0, 1], method='brentq')
        p = sol.root
    except ValueError:
        print("Error: Could not find a root in the given interval [0, 1].")
        return

    # Calculate the required value: floor(10000 * (1-p))
    one_minus_p = 1 - p
    value_to_floor = 10000 * one_minus_p
    final_answer = math.floor(value_to_floor)

    # Output the steps as requested
    print("The game has N races and 3N players.")
    print(f"For N = {N}, the number of players is {3*N}.")
    print("\nThe equilibrium probability 'p' is the root of the equation: p = 1 - (1 - p/N)^(3*N)")
    print("Substituting the values for this problem:")
    print(f"p = 1 - (1 - p/{N})^({3*N})")
    
    print(f"\nSolving this equation numerically yields p = {p:.10f}")
    
    print(f"\nNow we calculate the expression floor(10000 * (1-p)):")
    print(f"1 - p = {one_minus_p:.10f}")
    print(f"10000 * (1 - p) = {value_to_floor:.10f}")
    print(f"The floor of this value is {final_answer}")
    
    print(f"\nThe final answer is: {final_answer}")

# Execute the function
solve_for_p_and_calculate()