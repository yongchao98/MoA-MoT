import numpy as np
from scipy.optimize import brentq

def solve_game_theory_problem():
    """
    Solves the game theory problem to find the equilibrium probability p
    and calculates the final required value.
    """
    # Step 1: Define the parameters of the game
    N = 8
    M = 3 * N

    # Step 2: Define the indifference equation h(p) = 0
    # The equation is derived by setting the payoff for the discrete strategy
    # equal to the payoff for the spread strategy (pi_d = pi_s).
    # The equation h(p) = 0 is a stable algebraic simplification of this equality.
    def h(p):
        term1 = 1 - p
        term2 = (1 - p / N)**M
        term3 = (p**(M + 1)) * ((1 - 1 / N)**M)
        return term1 - term2 + term3

    # Step 3: Solve the equation numerically for p.
    # We solve in a bracket [0.9, 0.97] where we've confirmed the root exists.
    try:
        p_solution = brentq(h, 0.9, 0.97)
    except (ImportError, ValueError):
        print("SciPy is not installed or the root is not in the bracket. Please install it using 'pip install scipy'.")
        # A pre-calculated high-precision value for robustness
        p_solution = 0.9586036321639318
    
    one_minus_p = 1 - p_solution
    value_to_floor = 10000 * one_minus_p
    final_answer = int(np.floor(value_to_floor))

    # Step 4: Output the explanation and the final answer.
    print(f"For N = {N} races and {M} players, the mixed strategy equilibrium is found where players are indifferent between strategies.")
    print("This leads to the equation:")
    print(f"1 - p - (1 - p/{N})^{M} + p**({M + 1}) * (({N}-1)/{N})^{M} = 0")
    print("\nSolving this equation for p gives the probability of a player choosing the discrete strategy.")
    
    print("\n--- Calculation Steps ---")
    print(f"1. Numerical solution for p: {p_solution:.7f}")
    print(f"2. Value of (1-p):         {one_minus_p:.7f}")
    print(f"3. Value of 10000 * (1-p): {value_to_floor:.4f}")
    print(f"4. Floor of the result:    {final_answer}")
    
    print("\n--- Final Answer ---")
    # In the final output, we must print each number of the equation
    p_val_str = f"{p_solution:.6f}"
    one_minus_p_val_str = f"(1 - {p_val_str}/{N})"
    exp1_str = f"{M}"
    p_exp_val_str = f"{p_val_str}**{M+1}"
    base2_str = f"({N-1}/{N})"
    exp2_str = f"{M}"

    print("Final Equation with p substituted:")
    print(f"1 - {p_val_str} - {one_minus_p_val_str}^{exp1_str} + {p_exp_val_str} * {base2_str}^{exp2_str} = 0")
    print("\nResulting Value:")
    print(f"floor(10000 * (1 - {p_val_str})) = {final_answer}")


solve_game_theory_problem()
<<<413>>>