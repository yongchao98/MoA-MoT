import numpy as np
from scipy.optimize import brentq

def solve_and_display():
    """
    This function solves for the mixed strategy equilibrium probability p in the described game,
    and then calculates the final requested value. It also prints the details of the
    equilibrium equation for verification.
    """
    # Game parameters
    N = 8
    M = 3 * N  # Total number of players

    def equilibrium_equation(p):
        """
        This function represents the equilibrium condition pi_D(p) - pi_S(p) = 0.
        A root of this function in the interval (0, 1) is the equilibrium probability.
        
        Args:
            p: The probability of a player choosing the discrete strategy.
        
        Returns:
            The difference between the payoffs of the discrete and sprinkle strategies.
        """
        # The solver operates on a continuous interval, but p must be in (0, 1).
        if p <= 1e-9 or p >= 1 - 1e-9:
            return np.sign(0.5 - p) * 1e9

        # Payoff for the Discrete Strategy (pi_D)
        # This is E[1/(K+1)] where K ~ Bin(M-1, p/N).
        # The simplified formula is (1 - (1 - p/N)^M) / (3*p).
        pi_D = (1 - (1 - p / N)**M) / (3 * p)

        # Payoff for the Sprinkle Strategy (pi_S)
        # This is N * P(no discrete players) * P(win tie-break vs other sprinklers).
        prob_no_discrete = (1 - p / N)**(M - 1)
        
        # P(win tie-break) is E[1/(L+1)] where L ~ Bin(M-1, 1-p).
        # The simplified formula is (1 - p^M) / (M * (1-p)).
        expected_win_vs_sprinkle = (1 - p**M) / (M * (1 - p))
        
        pi_S = N * prob_no_discrete * expected_win_vs_sprinkle
        
        return pi_D - pi_S

    # Find the root of the equilibrium equation in the interval (0, 1).
    # We use brentq for robust root-finding in a bracketed interval.
    try:
        p_equilibrium = brentq(equilibrium_equation, 1e-9, 1 - 1e-9)
    except ValueError:
        print("Could not find a solution in the interval (0, 1).")
        return

    # Display the results as requested
    print("The mixed strategy equilibrium is found by solving pi_D(p) = pi_S(p) for p.")
    print(f"For N = {N}, the number of players is {M}.")
    print(f"The solved equilibrium probability is p = {p_equilibrium:.6f}\n")

    # Verify the equation with the solved value of p
    p = p_equilibrium
    
    # Calculate LHS (pi_D)
    lhs_val = (1 - (1 - p / N)**M) / (3 * p)
    
    # Calculate RHS (pi_S)
    prob_no_discrete_val = (1 - p / N)**(M - 1)
    expected_win_vs_sprinkle_val = (1 - p**M) / (M * (1 - p))
    rhs_val = N * prob_no_discrete_val * expected_win_vs_sprinkle_val

    print("Verification of the equilibrium equation with the solved p:")
    print(f"LHS = (1 - (1 - {p:.6f} / {N})**{M}) / (3 * {p:.6f}) = {lhs_val:.6f}")
    print(f"RHS = {N} * (1 - {p:.6f} / {N})**{M-1} * (1 - {p:.6f}**{M}) / ({M} * (1 - {p:.6f})) = {rhs_val:.6f}\n")

    # Calculate the final answer
    final_result = np.floor(10000 * (1 - p_equilibrium))
    
    print(f"The value of (1-p) is: {1 - p_equilibrium:.6f}")
    print(f"The value of 10000 * (1-p) is: {10000 * (1 - p_equilibrium):.6f}")
    print(f"The final result is floor(10000 * (1-p)) = {int(final_result)}")

solve_and_display()