import math

def solve_and_calculate_ev():
    """
    This function solves the problem by following the plan outlined above.
    """
    # Step 1: Determine the value of h.
    # The probability of traversing f, p_f, follows p_f = x + (1-x) * p_f^2.
    # Rearranging gives (1-x)p_f^2 - p_f + x = 0.
    # The solutions for p_f are 1 and x/(1-x).
    # The probability of success converges to 1 if and only if x > 0.5.
    # "h" is the smallest value of x where success is guaranteed, which is the boundary case.
    h = 0.5

    # Step 2: Determine the probabilities and success rate for g.
    # The three paths in g have probabilities p1, p2, p3.
    # p1 = h = 0.5 (direct success)
    # p2 = probability of hole (failure)
    # p3 = probability of chain of six g's
    # We are given p1 + p2 + p3 = 1 and p2 = 2 * p3.
    # 0.5 + 2*p3 + p3 = 1  =>  3*p3 = 0.5  =>  p3 = 1/6.
    p_chain_g = 1/6

    # The success probability of g, p_g, follows the equation:
    # p_g = (prob_direct * success_rate) + (prob_hole * success_rate) + (prob_chain * success_rate)
    # p_g = (h * 1) + (p_hole * 0) + (p_chain_g * p_g^6)
    # p_g = 0.5 + (1/6) * p_g^6
    # We solve this numerically for p_g using a few iterations.
    p_g = 0.5  # Initial guess
    for _ in range(10):  # 10 iterations is sufficient for convergence
        p_g = 0.5 + (1/6) * (p_g ** 6)

    # Step 3: Determine the success probability of k.
    # k is a chain of four instances of g.
    p_k = p_g ** 4

    # Step 4: Calculate the game's expected value.
    # The number of successes S in 100 traversals of k follows a binomial distribution
    # B(n, p) with n=100 and p=p_k.
    n = 100
    p = p_k

    # The opponent wins if S < 6. This is our probability of losing.
    # P(lose) = P(S < 6) = P(S=0) + P(S=1) + ... + P(S=5)
    p_lose = 0
    for i in range(6):
        # Binomial probability formula: C(n, i) * p^i * (1-p)^(n-i)
        term = math.comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
        p_lose += term
    
    p_win = 1.0 - p_lose

    # Expected Value = (value_win * p_win) + (value_lose * p_lose)
    # EV = ($1 * p_win) + (-$1 * p_lose)
    expected_value = (1 * p_win) + (-1 * p_lose)

    # Output the results as requested.
    print(f"The probability of successfully traversing a single 'k' structure is: {p:.6f}")
    print(f"The probability of you winning the bet (S >= 6) is: {p_win:.6f}")
    print(f"The probability of you losing the bet (S < 6) is: {p_lose:.6f}")
    print("\nThe final equation for the expected value is:")
    print(f"EV = ($1 * {p_win:.6f}) + (-$1 * {p_lose:.6f})")
    print(f"\nYour expected value of playing one game is ${expected_value:.2f}")

solve_and_calculate_ev()