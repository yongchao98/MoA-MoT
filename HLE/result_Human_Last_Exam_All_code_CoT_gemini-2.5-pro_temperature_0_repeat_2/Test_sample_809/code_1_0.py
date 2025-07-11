import scipy.optimize
import scipy.stats

def solve_and_calculate():
    """
    This function follows the plan to solve the problem step-by-step.
    1. Determines h.
    2. Solves for P(g).
    3. Calculates P(k).
    4. Computes the expected value of the bet.
    """
    # Step 1: Determine h
    # The probability of traversing f, P(f), is given by the smaller root of
    # P(f) = x + (1-x)P(f)^2. The roots are 1 and x/(1-x).
    # So, P(f) = min(1, x/(1-x)).
    # For traversal to be guaranteed, P(f) must be 1. This requires x/(1-x) >= 1,
    # which means x >= 0.5. The smallest such x is h.
    h = 0.5

    # Step 2: Determine P(g)
    # The probabilities for the paths in g are:
    # P_direct = h = 0.5
    # P_chain = 1/6
    # P_hole = 1/3
    # The probability of traversing g, P(g), is given by the recursive formula:
    # P(g) = 0.5 * 1 (direct path) + (1/3) * 0 (hole) + (1/6) * P(g)^6 (chain path)
    # This simplifies to P(g) = 0.5 + (1/6)P(g)^6.
    # We need to find the root of the equation f(p) = p^6 - 6p + 3 = 0 in [0, 1].
    def g_equation(p):
        return p**6 - 6*p + 3

    # Find the root using a numerical solver.
    p_g_solution = scipy.optimize.root_scalar(g_equation, bracket=[0, 1], method='brentq')
    p_g = p_g_solution.root

    # Step 3: Determine P(k)
    # k is a chain of 4 instances of g.
    p_k = p_g**4

    # Step 4: Calculate the expected value of the bet
    # The number of successes S in n=100 trials follows a binomial distribution B(n, p).
    n = 100
    p_success = p_k

    # The opponent wins if S < 6, which is S <= 5.
    # Your expected value is EV = P(S >= 6) - P(S <= 5) = 1 - 2 * P(S <= 5).
    # We calculate P(S <= 5) using the binomial CDF.
    num_successes_threshold = 5
    prob_opponent_wins = scipy.stats.binom.cdf(num_successes_threshold, n, p_success)

    # Final calculation
    expected_value = 1 - 2 * prob_opponent_wins

    # Print the results as requested
    print(f"The probability of the direct path in g is h = {h}")
    print(f"The probability of traversing a single g is P(g) = {p_g:.6f}")
    print(f"The probability of traversing k is P(k) = P(g)^4 = {p_k:.6f}")
    print(f"The number of trials is n = {n}")
    print(f"The opponent wins if successes are less than {num_successes_threshold + 1}")
    print("\nCalculating the expected value:")
    print(f"EV = 1 - 2 * P(Successes <= {num_successes_threshold})")
    print(f"The final equation is: EV = 1 - 2 * {prob_opponent_wins:.6f}")
    final_ev = 1 - 2 * prob_opponent_wins
    print(f"The result is: EV = {final_ev:.6f}")
    print(f"\nYour expected value of playing one game is ${final_ev:.2f}.")

solve_and_calculate()