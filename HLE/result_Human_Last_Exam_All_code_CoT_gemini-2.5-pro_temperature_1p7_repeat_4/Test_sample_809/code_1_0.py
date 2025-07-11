import math
import scipy.optimize
import scipy.stats

def solve_and_calculate():
    """
    This function carries out the step-by-step calculation to find the expected value.
    """
    # Part 1: Determine h
    # The smallest value of x for which traversing f is guaranteed (P(f)=1) is 0.5.
    h = 0.5

    # Part 2: Determine the probability of successfully traversing g, p_g
    # Probabilities of paths in g: p_direct=0.5, p_chain=1/6, p_hole=1/3.
    # The success probability p_g is defined by: p_g = 0.5 + (1/6) * p_g^6.
    # We solve the equivalent polynomial equation: p_g^6 - 6*p_g + 3 = 0.
    def g_equation(p):
        return p**6 - 6 * p + 3

    # Find the root between 0.5 and 1.
    p_g = scipy.optimize.brentq(g_equation, 0.5, 1)

    # Part 3: Determine the probability of successfully traversing k, p_k
    # k is a chain of four g's.
    p_k = p_g**4

    # Part 4: Calculate the probability of the opponent winning
    # n=100 trials, opponent wins if successes S < 6 (S <= 5).
    n_trials = 100
    k_successes_for_opponent_win = 5
    p_opponent_wins = scipy.stats.binom.cdf(k_successes_for_opponent_win, n_trials, p_k)

    # Part 5: Calculate your expected value
    # EV = 1 * P(You Win) - 1 * P(Opponent Wins) = 1 - 2 * P(Opponent Wins)
    expected_value = 1 - 2 * p_opponent_wins

    # Round the final result to the nearest cent.
    final_ev = round(expected_value, 2)
    
    print("The final calculation is based on the following equation and values:")
    print(f"Expected Value = 1 - 2 * P(Successes < 6)")
    print(f"P(Successes < 6) is the Binomial CDF for k={k_successes_for_opponent_win}, n={n_trials}, p=p_k")
    print(f"where p_k = p_g^4")
    print(f"and p_g is the solution to p_g = 0.5 + (1/6)*p_g^6\n")

    print(f"Calculated p_g = {p_g:.6f}")
    print(f"Calculated p_k = {p_k:.6f}")
    print(f"Calculated P(Opponent Wins) = {p_opponent_wins:.6f}\n")
    print(f"Final Equation:")
    print(f"Expected Value = 1 - 2 * {p_opponent_wins:.6f}")
    print(f"Expected Value = {expected_value:.6f}")
    print(f"\nYour expected value of playing one game is ${final_ev:.2f}")

solve_and_calculate()