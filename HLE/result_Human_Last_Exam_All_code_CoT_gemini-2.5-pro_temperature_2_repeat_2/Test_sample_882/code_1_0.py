import math
from scipy.optimize import brentq

def solve_equilibrium_p():
    """
    Solves for the mixed strategy equilibrium probability p and computes the final result.
    """
    N = 8
    M = 3 * N

    def payoff_discrete(p):
        """
        Calculates the expected payoff for a player using the discrete strategy
        against a population playing the mixed strategy with probability p.
        """
        # A player using the discrete strategy competes only with other discrete players.
        # This formula calculates the expected payoff, including tie-breaks.
        # Use a limit for p=0 to avoid division by zero.
        if p < 1e-9:
            return 1.0
        return (N / (M * p)) * (1 - (1 - p / N)**M)

    def payoff_split2(p):
        """
        Calculates the expected payoff for a player using the split-2 strategy.
        This payoff is modeled as the probability of entering at least one race
        that is not chosen by any discrete-strategy players.
        """
        num_opponents = M - 1
        return 2 * (1 - p / N)**num_opponents - (1 - 2 * p / N)**num_opponents

    def difference_function(p):
        """
        The function whose root we want to find for the equilibrium condition.
        """
        return payoff_discrete(p) - payoff_split2(p)

    # Analysis shows the equilibrium p is between 0.5 and 1.
    # We use the Brent's method to find the root of the difference function.
    try:
        p_equilibrium = brentq(difference_function, 0.5, 1.0)
    except ValueError:
        print("Error: Could not find the equilibrium point in the given interval.")
        return

    # Calculate the final value as requested by the problem.
    result = math.floor(10000 * (1 - p_equilibrium))
    
    print(result)

if __name__ == "__main__":
    solve_equilibrium_p()