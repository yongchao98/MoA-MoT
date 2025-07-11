import math

def solve_and_calculate_ev():
    """
    This function encapsulates the entire logic to solve the problem.
    1. Finds the probability of traversing g, P(g).
    2. Calculates the probability of traversing k, P(k).
    3. Models the bet as a binomial distribution.
    4. Calculates the expected value of the bet.
    """

    # Step 1: Determine h.
    # The probability of traversing f is P(f) = x + (1-x)P(f)^2.
    # The roots for P(f) are 1 and x/(1-x).
    # For P(f) to be guaranteed to be 1, we need x/(1-x) >= 1, which means x >= 0.5.
    # The smallest such value is h = 0.5.
    h = 0.5

    # Step 2: Find P(g) by solving p = 0.5 + (1/6)p^6.
    # This is equivalent to finding the root of f(p) = p^6 - 6p + 3 = 0.
    # We use a numerical bisection method for high precision.
    def g_equation(p):
        return p**6 - 6 * p + 3

    low = 0.5
    high = 1.0
    # 100 iterations are sufficient for high precision.
    for _ in range(100):
        mid = (low + high) / 2
        # The function is decreasing in the interval [0,1]
        if g_equation(mid) > 0:
            low = mid
        else:
            high = mid
    p_g = (low + high) / 2

    # Step 3: Calculate P(k), the success probability for one trial.
    # k is a chain of four g's.
    p_k = p_g**4

    # Step 4: Calculate the probability of losing the bet, P(X < 6).
    # This is a binomial distribution with n=100 trials and p=p_k success probability.
    n = 100
    # The opponent wins if X < 6, so we lose.
    prob_lose = 0
    for k in range(6):  # Sum probabilities for k = 0, 1, 2, 3, 4, 5
        # Binomial PMF: C(n, k) * p^k * (1-p)^(n-k)
        term = math.comb(n, k) * (p_k**k) * ((1 - p_k)**(n - k))
        prob_lose += term

    # Step 5: Calculate the final expected value.
    # E = (P(win) * $1) + (P(lose) * -$1)
    # P(win) = 1 - prob_lose
    # E = (1 - prob_lose) - prob_lose = 1 - 2 * prob_lose
    prob_win = 1 - prob_lose
    expected_value = prob_win - prob_lose

    # Print the final output as requested.
    print(f"The probability of you winning (X >= 6) is: {prob_win:.6f}")
    print(f"The probability of you losing (X < 6) is: {prob_lose:.6f}")
    print(f"Your expected value is ({prob_win:.6f}) * $1 + ({prob_lose:.6f}) * (-$1) = ${expected_value:.2f}")

solve_and_calculate_ev()