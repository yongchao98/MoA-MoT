import math

def solve_probability():
    """
    This function calculates the expected value of the game as described.
    """

    # Part 1: Find h
    # The probability of traversing f is P(f) = x + (1-x)P(f)^2.
    # For a guaranteed traversal, P(f) must be 1. This is possible if x >= 0.5.
    # h is the smallest such value.
    h = 0.5

    # Part 2: Find P(g), the probability of traversing g.
    # The probabilities for g's paths are:
    # p1 (direct success) = h = 0.5
    # p2 (hole) = 1/3
    # p3 (recursion) = 1/6
    # The traversal probability p_g satisfies: p_g = 0.5 + (1/6) * p_g^6
    # This is equivalent to finding the root of y^6 - 6y + 3 = 0 in [0, 1].
    
    def g_equation(y):
        return y**6 - 6*y + 3

    # Use bisection method to find the root p_g
    low, high = 0.0, 1.0
    for _ in range(100):  # 100 iterations for high precision
        mid = (low + high) / 2
        if g_equation(mid) > 0:
            low = mid
        else:
            high = mid
    p_g = (low + high) / 2

    # Part 3: Find P(k), the probability of traversing k.
    # k is a chain of 4 instances of g.
    p_k = p_g**4

    # Part 4: Calculate the probability of the opponent winning.
    # This is a binomial distribution with n=100 trials and success probability p=p_k.
    # The opponent wins if the number of successes is less than 6.
    n_trials = 100
    p_success = p_k
    
    prob_opponent_wins = 0.0
    for i in range(6):  # Sum probabilities for i = 0, 1, 2, 3, 4, 5 successes
        try:
            # Binomial probability: C(n, i) * p^i * (1-p)^(n-i)
            term = math.comb(n_trials, i) * (p_success**i) * ((1 - p_success)**(n_trials - i))
            prob_opponent_wins += term
        except ValueError:
            # This can happen if intermediate values are too large, but unlikely here.
            continue

    # Part 5: Calculate the expected value of the game.
    # E = P(win) * $1 + P(lose) * (-$1) = (1 - P(lose)) - P(lose) = 1 - 2 * P(lose)
    expected_value = 1 - 2 * prob_opponent_wins

    # Final Output: Print the final equation with calculated numbers and the result.
    print(f"The expected value is calculated by the equation: 1 - 2 * P(opponent wins)")
    print(f"1 - 2 * {prob_opponent_wins:.6f} = {expected_value:.6f}")
    print(f"Your expected value of playing one game is ${expected_value:.2f}.")

solve_probability()