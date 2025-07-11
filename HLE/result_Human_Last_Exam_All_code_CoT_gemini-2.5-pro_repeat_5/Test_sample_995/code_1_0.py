import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

def solve():
    """
    This function executes the plan to solve the problem numerically.
    """

    # 1. Define the probability of A winning given r and F, assuming A's optimal strategy.
    def p_a_wins_given_r(r, F):
        """
        Calculates the probability of agent A winning given the target's radial distance r
        and agent B's move distance F.
        """
        # If r < F/2, A can choose a move D -> 0 and win with probability 1.
        if r < F / 2.0:
            return 1.0
        else:
            # If r >= F/2, A's optimal move D gives a win probability based on arccos.
            # We must handle potential floating point inaccuracies at the boundaries.
            arg_sqrt = F * (2 * r - F)
            if arg_sqrt < 0:
                arg_sqrt = 0
            
            # C is the argument of arccos
            C = np.sqrt(arg_sqrt) / r
            if C > 1.0:
                C = 1.0
            
            return np.arccos(C) / np.pi

    # 2. Define the total probability of A winning for a given F by integrating over r.
    def p_a_wins(F):
        """
        Calculates the total probability of agent A winning for a given F.
        The probability density function for r is f(r) = 2r.
        P(F) = Integral[ P(A wins|r, F) * 2r dr ] from 0 to 1.
        """
        if F <= 0:
            return 0.5 # A known special case
        if F >= 2: # B moves outside the circle
            return 1.0
            
        # The integral is split at r = F/2 based on A's strategy.
        # Part 1: Integral from 0 to F/2. The integrand is 1.0 * 2r.
        prob_part1 = (F / 2.0)**2

        # Part 2: Integral from F/2 to 1. The integrand is p_a_wins_given_r(r, F) * 2r.
        integrand = lambda r: p_a_wins_given_r(r, F) * 2 * r
        prob_part2, _ = quad(integrand, F / 2.0, 1.0)
        
        return prob_part1 + prob_part2

    # 3. Use a numerical optimizer to find F in [0, 1] that minimizes P(A wins).
    result = minimize_scalar(p_a_wins, bounds=(0, 1), method='bounded')

    # 4. Extract the minimized probability.
    min_prob_A_wins = result.fun

    # 5. Calculate the final answer as floor(1 / P_min).
    final_value = np.floor(1 / min_prob_A_wins)

    # 6. Print the components of the final calculation and the result.
    # The final equation is: Result = floor(1 / P_min)
    print(f"The optimal value of F is approximately: {result.x}")
    print(f"The minimized probability P(A wins) is: {min_prob_A_wins}")
    print(f"The value of 1 / P(A wins) is: {1 / min_prob_A_wins}")
    print(f"The final result, floor(1 / P(A wins)), is: {final_value}")

solve()