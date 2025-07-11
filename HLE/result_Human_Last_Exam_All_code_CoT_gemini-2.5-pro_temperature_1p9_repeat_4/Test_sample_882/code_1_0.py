import math
import scipy.optimize

def solve_game_equilibrium():
    """
    Solves for the mixed strategy equilibrium probability 'p' in the described game.
    """
    N = 8
    PLAYERS = 3 * N
    OPPONENTS = PLAYERS - 1

    def pi_P(p):
        """
        Calculates the expected payoff for a player using the Discrete Strategy.
        """
        # Using a small epsilon to avoid division by zero or domain errors at p=0.
        if p <= 1e-9:
            # The limit as p -> 0 is 1, by L'Hopital's rule.
            return 1.0
        # The formula for the payoff is derived as (1 - (1 - p/N)**PLAYERS) / (3 * p)
        return (1.0 - (1.0 - p / N)**PLAYERS) / (3.0 * p)

    def pi_M(p):
        """
        Calculates the expected payoff for a player using the Spreading Strategy.
        """
        total = 0.0
        if p <= 1e-9:
            # The limit as p -> 0 is 1.
            return 1.0

        for k in range(1, N + 1):
            base = 1.0 - p * k / N
            # The base must be non-negative for real-valued powers, which holds for p in [0,1].
            term = ((-1.0)**(k - 1) * math.comb(N, k) * (base)**OPPONENTS)
            total += term
        return total

    def diff_function(p):
        """
        The difference between the two payoff functions. We need to find the root of this function.
        """
        return pi_P(p) - pi_M(p)

    # We need to find a bracketing interval [a, b] where diff_function(a) and diff_function(b) have opposite signs.
    # Testing shows the root is between 0.9 and 1.0.
    # diff_function(0.9) is slightly positive, diff_function(1.0) is negative.
    lower_bound = 0.9
    upper_bound = 1.0
    
    try:
        p_solution = scipy.optimize.brentq(diff_function, lower_bound, upper_bound)
    except ValueError:
        print("Error: Could not find a root in the given interval. The problem might be ill-defined.")
        return

    # Calculate the final result as requested by the problem.
    result = math.floor(10000 * (1.0 - p_solution))
    
    print("The probability p is: {:.6f}".format(p_solution))
    print("The value of 1 - p is: {:.6f}".format(1.0 - p_solution))
    print("The final result is: {}".format(result))
    
    # We must format the final output as requested in the problem.
    # No, the problem asks for the code to compute it. Printing the final value is fine.
    
solve_game_equilibrium()
