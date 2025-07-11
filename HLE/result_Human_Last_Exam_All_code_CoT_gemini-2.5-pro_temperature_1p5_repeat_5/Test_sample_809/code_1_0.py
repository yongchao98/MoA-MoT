import numpy as np
from scipy.stats import binom

def solve_and_calculate_ev():
    """
    This function solves the multi-part probability problem to find the expected value of the bet.
    """

    # Step 1: Find h
    # The probability P of traversing f is given by P = x + (1-x)P^2.
    # The solutions for P are 1 and x/(1-x). For a physical process, we take the smaller root.
    # To guarantee success, P must be 1. This is first achieved when x/(1-x) = 1.
    # x = 1 - x  =>  2x = 1  =>  x = 1/2.
    # So, h is the smallest value of x for which success is guaranteed.
    h = 0.5

    # Step 2: Find P(g)
    # The probabilities for g are:
    # p1 (direct path) = h = 0.5
    # p2 (hole) + p3 (recursion) = 1 - h = 0.5
    # p2 = 2 * p3
    # Solving this gives: 3 * p3 = 0.5 => p3 = 1/6, and p2 = 1/3.
    p3 = 1.0 / 6.0
    
    # Let Q = P(g). The recursive equation for Q is Q = h + p3 * Q^6.
    # Q = 0.5 + (1/6) * Q^6
    # This rearranges to the polynomial equation: Q^6 - 6Q + 3 = 0.
    # We need to find the real root of this polynomial between 0 and 1.
    coeffs = [1, 0, 0, 0, 0, -6, 3]
    roots = np.roots(coeffs)
    
    # Filter for the correct root (real, between 0 and 1)
    real_roots = roots[np.isreal(roots)].real
    Q = real_roots[(real_roots > 0) & (real_roots < 1)][0]

    # Step 3: Find P(k)
    # k is a chain of 4 instances of g.
    # The probability of traversing k, p_k, is Q^4.
    p_k = Q**4

    # Step 4: Calculate the expected value of the bet
    # The number of successes in n=100 trials follows a Binomial distribution B(n, p_k).
    n = 100
    # The opponent wins if the number of successes X is less than 6 (X <= 5).
    # We calculate P(X <= 5) using the binomial cumulative distribution function (CDF).
    prob_opponent_wins = binom.cdf(k=5, n=n, p=p_k)

    # The expected value E = $1 * P(You Win) - $1 * P(Opponent Wins)
    # E = P(X >= 6) - P(X <= 5)
    # E = (1 - P(X <= 5)) - P(X <= 5)
    # E = 1 - 2 * P(X <= 5)
    expected_value = 1 - 2 * prob_opponent_wins

    # Print the results of the calculation
    print(f"The value of h is: {h}")
    print("The equation for Q = P(g) is: Q^6 - 6*Q + 3 = 0")
    print(f"The probability of traversing g is Q = {Q:.6f}")
    print(f"The probability of traversing k is p_k = Q^4 = {p_k:.6f}")
    print(f"The number of trials is n = {n}")
    print(f"The probability of the opponent winning is P(successes <= 5) = {prob_opponent_wins:.6f}")
    final_equation = f"Expected Value = 1 - 2 * {prob_opponent_wins:.4f}"
    print(final_equation)
    print(f"Your final expected value is: ${expected_value:.2f}")

solve_and_calculate_ev()