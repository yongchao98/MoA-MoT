import numpy as np
from scipy.optimize import fsolve

def solve_tshirt_game_equilibrium():
    """
    Solves for the mixed-strategy Nash Equilibrium probability 'p' in the t-shirt game
    for N=8 and calculates the final requested value.
    """
    # Step 1: Define game parameters
    N = 8.0
    TOTAL_PLAYERS = 3.0 * N
    M = TOTAL_PLAYERS - 1.0  # Number of other players

    # Step 2: Define the constant C1 from the equilibrium equation
    C1 = (N - 1.0) * np.power((N - 1.0) / N, M)

    # Step 3: Define the function representing the equilibrium equation f(p) = 0
    # The equation is: pi_D(p) * (1 + C1) - ((1 - p/N)^M + C1) = 0
    def equilibrium_equation(p):
        # Handle the edge case of p=0, where pi_D requires L'Hopital's rule, giving pi_D(0)=1
        if p < 1e-9:
            p = 1e-9 # Use a small epsilon to avoid division by zero
        
        # Payoff for the discrete strategy
        pi_D = (1.0 - np.power(1.0 - p / N, M + 1.0)) / (3.0 * p)
        
        # Left-hand side of the rearranged equation
        lhs = pi_D * (1.0 + C1)
        
        # Right-hand side of the rearranged equation
        rhs = np.power(1.0 - p / N, M) + C1
        
        return lhs - rhs

    print("The equilibrium equation is derived from the condition: π_D(p) * (1 + C₁) = (1 - p/N)^M + C₁")
    print(f"For N={int(N)}, we have M={int(M)}.")
    # The term (1-p/N)^(M+1) is also expressible as (1-p/N)^(3N)
    print(f"π_D(p) = (1 - (1 - p/{N})^{{{int(M)+1}}}) / (3*p)")
    print(f"C₁ = ({int(N)}-1) * (1 - 1/{int(N)})^{{{int(M)}}} ≈ {C1:.6f}")
    
    # Step 4: Solve the equation numerically for p
    # We need a reasonable starting guess for the solver.
    initial_guess = 0.5
    p_solution = fsolve(equilibrium_equation, initial_guess)
    p = p_solution[0]
    
    print(f"\nSolved for p, the probability of playing the discrete strategy.")
    print(f"The value of p is approximately: {p:.6f}")

    # Step 5: Calculate and print the final result
    one_minus_p = 1 - p
    result = np.floor(10000 * one_minus_p)
    
    print(f"\nThe value of (1-p) is: {one_minus_p:.6f}")
    print(f"The value of 10000 * (1-p) is: {10000 * one_minus_p:.6f}")
    print(f"The final answer, floor(10000 * (1-p)), is:")
    print(int(result))

solve_tshirt_game_equilibrium()