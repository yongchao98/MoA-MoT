import math
import sys
import subprocess

# Ensure scipy is installed
try:
    from scipy.optimize import brentq
except ImportError:
    print("SciPy library not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scipy"])
    from scipy.optimize import brentq

def solve_game_theory_problem():
    """
    Solves the game theory problem to find the value of floor(10000 * (1-p)).
    """
    # Game parameters from the problem description
    N = 8
    THREE_N = 3 * N
    THREE_N_MINUS_1 = 3 * N - 1
    
    # For N=8, the optimal deviating strategy is to split fuel across m=4 races.
    OPTIMAL_M = 4

    print("The problem is set for N=8 races and 3N=24 players.")
    print(f"The equilibrium is a mix between a 'discrete' strategy and a 'split' strategy over m={OPTIMAL_M} races.")
    print("We need to solve for p where Payoff(discrete) = Payoff(split).")
    print("\n--- The Equation ---")
    
    # Define and print the formulas for the payoffs
    lhs_formula = f"(1 - (1 - p/{N})^{THREE_N}) / (3*p)"
    print(f"Payoff(discrete) = {lhs_formula}")
    
    rhs_formula_parts = []
    for j in range(1, OPTIMAL_M + 1):
        sign = "-" if (j - 1) % 2 else "+"
        if j == 1: sign = ""
        term = f"{sign} {math.comb(OPTIMAL_M, j)}*(1 - p*{j}/{N})^{THREE_N_MINUS_1}"
        rhs_formula_parts.append(term)
    
    print(f"Payoff(split, m={OPTIMAL_M}) = {' '.join(rhs_formula_parts)}")
    print("--------------------\n")

    def payoff_discrete(p):
        """
        Calculates the payoff for the discrete strategy.
        p: probability of playing discrete.
        """
        if p == 0:
            # The limit as p -> 0 is 1.
            return 1.0
        return (1 - (1 - p / N)**THREE_N) / (3 * p)

    def payoff_split(p, m):
        """
        Calculates the payoff for the split-m strategy.
        p: probability of playing discrete.
        m: number of races to split fuel into.
        """
        total = 0
        for j in range(1, m + 1):
            term = ((-1)**(j - 1) * math.comb(m, j) * (1 - p * j / N)**THREE_N_MINUS_1)
            total += term
        return total

    def indifference_equation(p):
        """
        The function g(p) = Payoff(split) - Payoff(discrete).
        We want to find the root p where g(p) = 0.
        """
        # A small value of p can cause precision issues, but the root is not near 0.
        if p <= 1e-12:
            return 1.0 # Function is positive for small p > 0

        p_d = payoff_discrete(p)
        p_s = payoff_split(p, OPTIMAL_M)
        return p_s - p_d

    # We solve for p in the interval (0, 1). Based on analysis, the root
    # is between 0.9 and 1.0. brentq is an efficient root-finding algorithm.
    try:
        p_solution = brentq(indifference_equation, 0.9, 1.0)
    except ValueError:
        print("Error: Could not find the root in the specified interval [0.9, 1.0].")
        print("This may indicate an issue with the model or interval.")
        return

    # Calculate the final result as requested by the problem
    one_minus_p = 1 - p_solution
    final_value = math.floor(10000 * one_minus_p)

    print(f"Numerically solving for p gives a value of approximately {p_solution:.6f}.")
    print(f"The value of 1 - p is {one_minus_p:.6f}.")
    print(f"The expression 10000 * (1-p) is approximately {10000 * one_minus_p:.4f}.")
    print(f"The final result, floor(10000 * (1-p)), is:")
    print(final_value)

solve_game_theory_problem()
<<<871>>>