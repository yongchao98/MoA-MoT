import math

def find_root_bisection(func, a, b, tol=1e-12):
    """
    Finds the root of a function using the bisection method.
    """
    if func(a) * func(b) >= 0:
        raise ValueError("Function has the same sign at interval ends.")
    while (b - a) / 2.0 > tol:
        midpoint = (a + b) / 2.0
        if func(midpoint) == 0:
            return midpoint
        elif func(a) * func(midpoint) < 0:
            b = midpoint
        else:
            a = midpoint
    return (a + b) / 2.0

def binomial_cdf(k_max, n, p):
    """
    Calculates the cumulative distribution function P(X <= k_max) for a binomial distribution.
    """
    cumulative_prob = 0.0
    for i in range(k_max + 1):
        # P(X=i) = C(n, i) * p^i * (1-p)^(n-i)
        term = math.comb(n, i) * (p**i) * ((1-p)**(n-i))
        cumulative_prob += term
    return cumulative_prob

def solve():
    """
    Solves the entire problem step-by-step.
    """
    # Part 1: Define h
    h = 0.5

    # Part 2: Define probabilities for g and find P(g)
    prob_g_direct = h
    # We have prob_g_hole = 2 * prob_g_chain and h + prob_g_hole + prob_g_chain = 1
    # 0.5 + 3 * prob_g_chain = 1 => prob_g_chain = 0.5 / 3
    prob_g_chain = 1/6
    
    # Equation for p_g is: p_g = prob_g_direct * 1 + prob_g_chain * p_g^6
    # p_g = 0.5 + (1/6) * p_g^6
    # Find root of: y - 0.5 - (1/6)*y^6 = 0. We know root is in (0.5, 1).
    def g_traversal_equation(y):
        return y - prob_g_direct - prob_g_chain * (y**6)

    p_g = find_root_bisection(g_traversal_equation, 0.5, 1.0)

    # Part 3: Define P(k)
    num_g_in_k = 4
    p_k = p_g**num_g_in_k

    # Part 4: Calculate expected value
    # The final equation for the expected value E is:
    # E = 1 * P(win) - 1 * P(lose) = 1 - 2 * P(lose)
    # P(lose) is the probability of having less than 6 successes in 100 trials.
    # P(lose) = P(X <= 5)
    num_trials = 100
    success_threshold = 5 # lose if successes <= 5
    
    p_lose = binomial_cdf(success_threshold, num_trials, p_k)
    
    expected_value = 1 - 2 * p_lose
    
    print("--- Step-by-step calculation ---")
    print(f"Value of h: {h}")
    print(f"Probabilities for g: Direct={prob_g_direct}, Chain={prob_g_chain:.4f}, Hole={2*prob_g_chain:.4f}")
    print(f"Traversal probability for g (p_g), root of y = {prob_g_direct} + ({prob_g_chain:.4f})*y^6, is: {p_g:.6f}")
    print(f"Traversal probability for k (p_k = p_g^{num_g_in_k}) is: {p_k:.6f}")
    print(f"Number of trials (n): {num_trials}")
    print(f"Opponent wins if successes < 6 (i.e., <= {success_threshold})")
    print(f"Your loss probability P(X <= {success_threshold}) is: {p_lose:.6f}")
    print("\nFinal expected value E = 1 - 2 * P(lose)")
    print(f"E = 1 - 2 * {p_lose:.6f}")
    print(f"E = {expected_value:.6f}")
    print(f"\nYour expected value is ${expected_value:.2f}")
    
solve()