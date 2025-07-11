import math

def solve():
    """
    This function calculates the expected value of the game as described.
    """
    # Step 1: Define h, the probability of the first path in g.
    # From the analysis of f, h is the smallest x for which success is guaranteed.
    # p = x + (1-x)p^2. Guaranteed success (p=1) happens for x >= 0.5.
    h = 0.5

    # Step 2: Define P(g), the probability of successfully traversing g.
    # The equation for q = P(g) is q = h + (1/6)q^6.
    # We will approximate q as h, so P(g) = 0.5.
    P_g = 0.5
    
    # Step 3: Define P(k), the probability of successfully traversing k.
    # k is a chain of 4 instances of g.
    # The probability is P_g to the power of 4.
    num_g_in_k = 4
    P_k = P_g ** num_g_in_k

    # Step 4: Calculate the expected value of the bet.
    # The number of trials for traversing k.
    n_trials = 100
    
    # The opponent bets on the number of successes being less than 6.
    # This means the opponent wins if successes are 0, 1, 2, 3, 4, or 5.
    win_threshold = 6
    
    # Calculate the probability of the opponent winning (your P_lose).
    # This is the cumulative distribution function (CDF) of the binomial distribution.
    # P(X <= 5) = sum from i=0 to 5 of P(X=i)
    p_lose = 0
    for i in range(win_threshold):
        # Binomial probability: C(n, i) * p^i * (1-p)^(n-i)
        prob_i = math.comb(n_trials, i) * (P_k**i) * ((1 - P_k)**(n_trials - i))
        p_lose += prob_i

    # Expected Value (EV) = 1 * P(win) - 1 * P(lose) = 1 - 2 * P(lose)
    expected_value = 1 - 2 * p_lose

    # The final equation is constructed by showing all the numerical parts.
    # We will show the calculation for P(X<=5) and the final EV calculation.
    print(f"h = 1/2")
    print(f"Probabilities for g are path 1 (direct success): {h}, path 2 (hole): {1/3:.3f}, path 3 (recursive): {1/6:.3f}")
    print(f"Probability of success for g, P(g), is approximated as {P_g}")
    print(f"Probability of success for k, P(k) = P(g)^4 = {P_g}^4 = {P_k}")
    print(f"The probability of your opponent winning is P(Successes < {win_threshold}), which is P(Successes <= {win_threshold-1})")
    print(f"P(Successes <= {win_threshold-1}) = {p_lose:.6f}")
    print(f"Your expected value = $1 * P(win) - $1 * P(lose) = 1 - 2 * P(lose)")
    print(f"Expected Value = 1 - 2 * {p_lose:.6f} = {expected_value:.6f}")
    print(f"Your expected value rounded to the nearest cent is ${expected_value:.2f}")

solve()