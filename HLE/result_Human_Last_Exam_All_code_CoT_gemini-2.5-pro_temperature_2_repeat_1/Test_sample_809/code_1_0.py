import math

def solve():
    """
    This function solves the multi-step probability problem to find the expected value of the bet.
    """

    # Step 1: Determine h
    # The traversal probability P(f) for a given probability x follows P(f) = x + (1-x)*P(f)^2.
    # We can rewrite this as a quadratic equation: (1-x)P(f)^2 - P(f) + x = 0.
    # The traversal is 'guaranteed' (P(f)=1) for the smallest value of x when the quadratic
    # has exactly one real root, which happens when the discriminant is 0.
    # Discriminant = (-1)^2 - 4*(1-x)*x = 1 - 4x + 4x^2 = (1-2x)^2.
    # (1-2x)^2 = 0  => x = 0.5.
    # Therefore, h, the smallest value of x for a guaranteed traversal, is 0.5.
    h = 0.5

    # Step 2: Determine the traversal probability of g, P(g)
    # The probabilities of the paths in g are:
    # p1 (direct path) = h = 0.5
    # The other two paths have probabilities p2 and p3, where p1+p2+p3=1 and p2=2*p3.
    # 0.5 + 2*p3 + p3 = 1  => 3*p3 = 0.5  => p3 = 1/6.
    # p2 (hole) = 2/6 = 1/3.
    # The traversal probability P(g) is given by the recursive formula:
    # P(g) = (p1 * 1) + (p2 * 0) + (p3 * P(g)^6) = 0.5 + (1/6) * P(g)^6.
    # We solve for P(g) using numerical iteration.
    p_g = 0.5  # Start with an initial guess
    for _ in range(10):  # 10 iterations are more than enough for convergence
        p_g = 0.5 + (1/6) * (p_g**6)

    # Step 3: Determine the traversal probability of k, P(k)
    # k is a chain of 4 instances of g.
    p_k = p_g**4

    # Step 4: Calculate the probability of losing the bet
    # This is a binomial distribution problem with n=100 trials and success probability p=p_k.
    n = 100
    p = p_k
    # The opponent wins if the number of successes is less than 6 (i.e., 0, 1, 2, 3, 4, or 5).
    # P(lose) = P(X < 6)
    p_lose = 0
    for i in range(6):
        # Binomial probability: C(n, i) * p^i * (1-p)^(n-i)
        p_lose += math.comb(n, i) * (p**i) * ((1 - p)**(n - i))

    # Step 5: Calculate the expected value of the game
    p_win = 1 - p_lose
    expected_value = (p_win * 1.0) + (p_lose * -1.0)

    # Print the results
    print(f"Step 1: The value of h is {h}")
    print(f"Step 2: The probability of traversing g, P(g), is approximately {p_g:.6f}")
    print(f"Step 3: The probability of traversing k, P(k), is approximately {p_k:.6f}")
    print("\n--- Bet Calculation ---")
    print(f"The probability of you winning the bet, P(X >= 6), is approximately {p_win:.6f}")
    print(f"The probability of you losing the bet, P(X < 6), is approximately {p_lose:.6f}")
    print("\n--- Final Expected Value ---")
    print("Expected Value = P(win) * $1 + P(lose) * -$1")
    print(f"Expected Value = {p_win:.6f} * $1 + {p_lose:.6f} * -$1 = ${expected_value:.2f}")

solve()