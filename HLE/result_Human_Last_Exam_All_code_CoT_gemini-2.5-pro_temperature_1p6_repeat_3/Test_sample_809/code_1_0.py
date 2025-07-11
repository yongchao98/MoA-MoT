import scipy.optimize
import scipy.stats

def solve_and_calculate():
    """
    This function solves the problem by following the plan laid out.
    1. Determines h.
    2. Solves for the probability of traversing g, P(g).
    3. Calculates the probability of traversing k, P(k).
    4. Calculates the expected value of the described game.
    """
    # Step 1: Find h
    # The probability of traversing f, P(f), satisfies the equation: P(f) = x + (1-x) * P(f)^2.
    # Rearranging, we get (1-x)P(f)^2 - P(f) + x = 0.
    # The solutions for P(f) are 1 and x/(1-x). For P(f) to be a valid probability <= 1,
    # the second solution is valid only when x <= 0.5.
    # The probability of success becomes guaranteed (P(f) = 1) when x >= 0.5.
    # The smallest value of x for which this is true is h.
    h = 0.5
    print("--- Step 1: Determine h ---")
    print("The probability of traversing system f is guaranteed to be 1 when x >= 0.5.")
    print(f"The smallest such value, h, is {h}.")

    # Step 2: Find the probability of traversing g, P(g)
    # For g, P(direct) = h = 0.5. The remaining probability is 1 - h = 0.5.
    # This is split between a path to a hole (P_hole) and a path to a chain of g's (P_chain).
    # We are given P_hole = 2 * P_chain.
    # So, 2*P_chain + P_chain = 0.5 => 3*P_chain = 0.5 => P_chain = 1/6.
    prob_g_chain = 1/6

    # The success probability P(g) must satisfy:
    # P(g) = (Prob of direct path) * 1 + (Prob of chain path) * P(g)^6
    # P(g) = 0.5 + (1/6) * P(g)^6
    # This is equivalent to finding the root of f(p) = p^6 - 6*p + 3 = 0.
    def g_equation(p):
        return p**6 - 6*p + 3

    # Solve for P(g) numerically in the interval [0, 1].
    solution = scipy.optimize.root_scalar(g_equation, bracket=[0, 1], method='brentq')
    p_g = solution.root
    print("\n--- Step 2: Determine the probability of traversing g, P(g) ---")
    print(f"The success probability of g, P(g), solves the equation: P(g) = {h} + {prob_g_chain:.6f} * P(g)^6")
    print(f"The numerical solution is P(g) = {p_g:.8f}")

    # Step 3: Find the probability of traversing k, P(k)
    # k is a chain of four instances of g.
    p_k = p_g**4
    print("\n--- Step 3: Determine the probability of traversing k, P(k) ---")
    print("k is a chain of 4 instances of g.")
    print(f"P(k) = P(g)^4 = ({p_g:.8f})^4 = {p_k:.8f}")

    # Step 4: Calculate the expected value of the game
    n_trials = 100
    p_success = p_k
    bet_threshold = 6

    # The opponent wins if the number of successes, S, is less than 6 (i.e., S <= 5).
    # We use the binomial cumulative distribution function (CDF).
    prob_opponent_wins = scipy.stats.binom.cdf(k=bet_threshold - 1, n=n_trials, p=p_success)
    prob_you_win = 1 - prob_opponent_wins

    # Expected Value E = ($1 * P(You Win)) - ($1 * P(Opponent Wins))
    expected_value = prob_you_win - prob_opponent_wins

    print("\n--- Step 4: Calculate the Expected Value of the Game ---")
    print(f"The game involves {n_trials} trials with success probability p = {p_success:.8f}.")
    print("Your opponent wins if successful traversals < 6.")
    print("You win if successful traversals >= 6.")
    
    print("\nThe final equation for your expected value is:")
    print(f"Expected Value = ($1 * P(You Win)) - ($1 * P(Opponent Wins))")
    print(f"Expected Value = ($1 * {prob_you_win:.8f}) - ($1 * {prob_opponent_wins:.8f})")
    print(f"Expected Value = ${expected_value:.8f}")

    print(f"\nRounded to the nearest cent, your expected value is: ${expected_value:.2f}")

solve_and_calculate()