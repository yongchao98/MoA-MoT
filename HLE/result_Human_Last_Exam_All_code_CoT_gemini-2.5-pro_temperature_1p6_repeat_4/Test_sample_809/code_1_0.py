import math

def solve_puzzle():
    """
    This script calculates the expected value of the game described.
    """
    
    # Step 1: Determine h.
    # The probability of traversing f is P(f) = x + (1-x)P(f)^2.
    # Success is guaranteed (P(f)=1) for the smallest x where this is the only valid outcome.
    # This occurs when x = 0.5.
    h = 0.5

    # Step 2: Determine probabilities for the paths in g and solve for P(g).
    # p1 = h = 0.5. Given p1+p2+p3=1 and p2=2*p3, we find p3=1/6.
    # The success probability of g, q=P(g), follows q = 0.5 + (1/6)q^6.
    # We solve the equation q^6 - 6q + 3 = 0 for q in [0,1].
    
    def g_equation(q):
        return q**6 - 6 * q + 3

    # We use a bisection method to find the root numerically.
    low, high = 0.5, 1.0
    for _ in range(100):  # 100 iterations provides high precision
        mid = (low + high) / 2
        if g_equation(mid) * g_equation(high) > 0:
            high = mid
        else:
            low = mid
    q = (low + high) / 2
    
    # Step 3: Calculate the success probability of k.
    # P(k) is the probability of traversing 4 g's in a chain.
    p_k = q**4

    # Step 4: Calculate the expected value of the bet.
    # This involves the binomial distribution B(n=100, p=p_k).
    # We lose if the number of successes Y is less than 6.
    n_trials = 100
    p_lose = 0.0
    
    # Sum the probabilities for Y = 0, 1, 2, 3, 4, 5.
    for k_successes in range(6):
        # Calculate the binomial probability: C(n, k) * p^k * (1-p)^(n-k)
        try:
            # math.comb is efficient and available in Python 3.8+
            comb = math.comb(n_trials, k_successes)
        except AttributeError:
            # Fallback for older Python versions
            comb = math.factorial(n_trials) / (math.factorial(k_successes) * math.factorial(n_trials - k_successes))
            
        term = comb * (p_k**k_successes) * ((1 - p_k)**(n_trials - k_successes))
        p_lose += term

    # Final EV equation: EV = 1 - 2 * P(lose)
    expected_value = 1 - 2 * p_lose

    # Output the intermediate values and the final result.
    print(f"The probability h is: {h}")
    print(f"The success probability for a single g traversal is q = {q}")
    print(f"The success probability for a single k traversal is P(k) = q^4 = {p_k}")
    print(f"The probability of you losing the bet is P(successes < 6) = {p_lose}")
    print(f"The expected value calculation is: 1 - 2 * {p_lose}")
    print(f"Your final expected value is: ${expected_value:.2f}")

solve_puzzle()