import math
from scipy.optimize import brentq

def solve_game_theory_problem():
    """
    Solves for the probability 'p' in the mixed strategy Nash Equilibrium
    and computes the final requested value.
    """
    # Game parameters
    N = 8
    TOTAL_PLAYERS = 3 * N
    
    # Constants for the equation
    M = 3 * N
    M_PLUS_1 = M + 1
    
    # The equation for p is derived from setting the expected utilities of the
    # 'discrete' and 'split' strategies equal: U_D(p) = U_S(p).
    # This simplifies to: 1 - p = (1 - p/N)^M - p^(M+1) * ((N-1)/N)^M
    
    # We define a function F(p) such that F(p) = 0 at the equilibrium.
    # F(p) = (1 - p/N)^M - p^(M+1) * ((N-1)/N)^M - (1 - p)
    c = ((N - 1) / N) ** M
    
    def F(p):
        term1 = (1 - p / N) ** M
        term2 = p ** M_PLUS_1 * c
        term3 = 1 - p
        return term1 - term2 - term3

    # Print the equation with the specific numbers for N=8
    print("The equilibrium condition leads to the following equation for p:")
    print(f"1 - p = (1 - p/{N})^{M} - p^{M_PLUS_1} * (({N}-1)/{N})^{M}")
    print(f"For N=8, this is: 1 - p = (1 - p/8)^24 - p^25 * (7/8)^24")
    print("-" * 20)

    # Solve for p numerically. We know p must be between 0 and 1.
    # Testing shows the root is close to 1.
    try:
        p_equilibrium = brentq(F, 0.9, 1.0)
    except (ImportError, ValueError):
        print("Scipy not found or root not in bracket. Please ensure scipy is installed.")
        # Fallback to a simpler method or exit if required. For this problem, we'll assume it exists.
        return

    # 'p' is the probability that a player in Nash equilibrium devotes all of their fuel to a single race.
    print(f"The calculated probability p is: {p_equilibrium:.6f}")

    # Calculate the final value as requested by the problem
    value = 10000 * (1 - p_equilibrium)
    final_answer = math.floor(value)

    print(f"\nThe value of 1 - p is: {1 - p_equilibrium:.6f}")
    print(f"The value of 10000 * (1 - p) is: {value:.4f}")
    print(f"The floor of 10000 * (1 - p) is: {final_answer}")
    
    # Final answer in the specified format
    print(f"\n<<<{final_answer}>>>")


solve_game_theory_problem()
