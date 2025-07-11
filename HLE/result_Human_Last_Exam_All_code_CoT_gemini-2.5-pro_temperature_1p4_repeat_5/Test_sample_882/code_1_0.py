import math
from scipy import optimize

def solve_game_theory_problem():
    """
    This function solves the problem as described.
    It finds the value of p for N=8 and then calculates the final result.
    """
    # The number of races, as given in the problem.
    N = 8
    
    # The total number of players.
    num_players = 3 * N

    # In a symmetric Nash Equilibrium, the expected payoff for any strategy
    # in the equilibrium support must be equal. The average payoff is N / (3N) = 1/3.
    # The expected payoff for the discrete strategy is E = (1 - (1 - p/N)^(3N)) / (3p).
    # Setting E = 1/3 gives the equation: 1 - (1 - p/N)^(3N) = p.
    # We define a function f(p) = 0 to solve for p.
    # f(p) = 1 - p - (1 - p/N)^(3N)
    
    # The numbers in the final equation are:
    constant_term = 1
    p_coefficient = -1
    power_base_denominator = N
    exponent = num_players
    
    # Define the function whose root we want to find.
    f = lambda p: 1 - p - (1 - p / power_base_denominator)**exponent
    
    # We look for a non-trivial solution for p in the interval (0, 1).
    # We use the Brent's method (brentq) for robust root-finding.
    # We search in a slightly reduced interval to avoid boundary issues.
    try:
        p_solution = optimize.brentq(f, 1e-9, 1 - 1e-9)
    except (ValueError, RuntimeError) as e:
        print(f"Could not find a solution for p: {e}")
        return

    # Calculate the components of the final answer.
    one_minus_p = 1 - p_solution
    value = 10000 * one_minus_p
    final_answer = math.floor(value)

    # Print the step-by-step derivation as per the plan.
    print(f"Step 1: The number of races is N = {N}, and players = {num_players}.")
    print("Step 2: The equilibrium condition leads to the equation for the probability p:")
    print(f"   1 - p - (1 - p/{power_base_denominator})^{exponent} = 0")
    print(f"Step 3: Solving this equation numerically for p.")
    print(f"   The probability p is: {p_solution:.6f}")
    print(f"Step 4: Calculating the final expression.")
    print(f"   1 - p = {one_minus_p:.6f}")
    print(f"   10000 * (1 - p) = {value:.6f}")
    print(f"   The floor of the result is: {final_answer}")

solve_game_theory_problem()