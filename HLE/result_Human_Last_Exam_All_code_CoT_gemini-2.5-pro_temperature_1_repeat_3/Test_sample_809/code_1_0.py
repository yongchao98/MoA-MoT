import math
from scipy.optimize import root_scalar

def solve_and_print_ev():
    """
    This function solves the problem by following the plan and prints the results.
    """
    # Step 1: Determine the value of h.
    # The probability of traversing f, P(f), follows the recursive equation:
    # P(f) = x + (1-x) * P(f)^2, which gives solutions P(f)=1 and P(f)=x/(1-x).
    # For success to be guaranteed (P(f)=1), we need x/(1-x) >= 1, so x >= 0.5.
    # The smallest value 'h' for which this holds is 0.5.
    h = 0.5

    # Step 2: Determine the success probability of g, P(g).
    # For g, Path 1 prob = h = 0.5. Let Path 3 prob = p3. Then Path 2 prob = 2*p3.
    # 0.5 + 2*p3 + p3 = 1 => 3*p3 = 0.5 => p3 = 1/6.
    # The recursive equation for P(g) is: p_g = 0.5 * 1 + (1/6) * p_g^6.
    # This requires finding the root of the equation: p^6 - 6p + 3 = 0.
    def g_equation(p):
        return p**6 - 6 * p + 3

    # We find the root numerically in the valid probability range (0, 1).
    # Testing shows the root is between 0.5 and 0.6.
    solution = root_scalar(g_equation, bracket=[0.5, 0.6])
    p_g = solution.root

    # Step 3: Determine the success probability of k, P(k).
    # k is a chain of four g's, so P(k) = P(g)^4.
    p_k = p_g**4

    # Step 4: Calculate the expected value of the bet.
    # This is a binomial distribution with n=100 trials and success probability p_k.
    n = 100
    # The opponent wins if successes < 6 (i.e., 0, 1, 2, 3, 4, or 5).
    prob_opponent_wins = 0.0
    for i in range(6):  # Sum probabilities for i = 0 to 5
        term_prob = math.comb(n, i) * (p_k**i) * ((1 - p_k)**(n - i))
        prob_opponent_wins += term_prob

    # Your win probability is the complement.
    prob_you_win = 1.0 - prob_opponent_wins

    # The expected value is E = ($1 * P(You Win)) + (-$1 * P(Opponent Wins)).
    expected_value = prob_you_win - prob_opponent_wins

    # Print all the numbers involved in the final calculation.
    print(f"The probability of the first path in g is h = {h}.")
    print(f"The success probability of traversing g, P(g), is the root of p^6 - 6p + 3 = 0, which is P(g) = {p_g:.8f}.")
    print(f"The success probability of traversing k is P(k) = P(g)^4 = {p_k:.8f}.")
    print(f"The number of trials is n = {n}.")
    print(f"The opponent wins if the number of successes is less than 6.")
    print(f"Your win probability is P(successes >= 6) = {prob_you_win:.8f}.")
    print(f"Your opponent's win probability is P(successes < 6) = {prob_opponent_wins:.8f}.")
    print("\nThe final equation for your expected value is:")
    print(f"Expected Value = ($1 * Your Win Probability) - ($1 * Opponent's Win Probability)")
    print(f"Expected Value = $1 * {prob_you_win:.8f} - $1 * {prob_opponent_wins:.8f}")
    print(f"Expected Value = ${expected_value:.8f}")
    print(f"\nRounded to the nearest cent, your expected value is: ${expected_value:.2f}")

solve_and_print_ev()