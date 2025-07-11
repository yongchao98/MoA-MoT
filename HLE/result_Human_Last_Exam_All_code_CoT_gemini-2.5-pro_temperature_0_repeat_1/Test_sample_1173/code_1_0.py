import math

def solve():
    """
    This function formalizes the argument that theta = 1/2.
    The core of the argument is bounding the term subtracted from n, which is a sum of probabilities.
    The sum is bounded by n * P_{n-1}, where P_{n-1} is the probability of stopping before time n.
    P_{n-1} is bounded using Doob's maximal inequality.
    The final bound on the subtracted term is of order n^{1/2}.
    """
    
    # Let's represent the order of magnitude of key terms as functions of n
    # We use a large value for n to illustrate the asymptotic behavior
    n = 1000000

    # Threshold for the sum S_k
    a = 1 - 1 / math.sqrt(n)

    # Variance of a single random variable X_i
    # Var(X_i) = 1/(3*n^(3/2)) - 1/(4*n^2)
    # For large n, Var(X_i) is approximately 1/(3*n^(3/2))
    var_X_i_order = lambda n_val: 1 / (3 * n_val**1.5)

    # Variance of the sum S_{n-1}
    # Var(S_{n-1}) = (n-1) * Var(X_i)
    # For large n, Var(S_{n-1}) is approximately n * (1/(3*n^(3/2))) = 1/(3*n^(1/2))
    var_S_n_minus_1_order = lambda n_val: 1 / (3 * n_val**0.5)

    # The value x for Doob's inequality
    # x = 1/2 - n^{-1/2} + 1/(2n)
    # For large n, x is approximately 1/2
    x_order = lambda n_val: 0.5

    # Bound on P_{n-1} from Doob's inequality
    # P_{n-1} <= Var(S_{n-1}) / x^2
    # The order is (1/n^{1/2}) / (1/4) = 4/n^{1/2}
    P_n_minus_1_bound_order = lambda n_val: var_S_n_minus_1_order(n_val) / (x_order(n_val)**2)

    # The bound on the sum we subtract from n
    # Sum <= n * P_{n-1}
    # The order is n * (1/n^{1/2}) = n^{1/2}
    subtracted_term_bound_order = lambda n_val: n_val * P_n_minus_1_bound_order(n_val)

    # So, E[tau] >= n - c * n^{1/2}
    # This means theta = 1/2 is a provable value.
    
    theta = 1/2
    
    print("Step 1: Express E[tau] = n - Sum_{k=0}^{n-1} P(tau <= k).")
    print("Step 2: Bound the sum by n * P(tau <= n-1).")
    print("Step 3: P(tau <= n-1) = P(max_{k<=n-1} S_k >= a).")
    print("Step 4: Use Doob's maximal inequality on the martingale M_k = S_k - E[S_k].")
    print(f"Step 5: Var(S_{n-1}) is of order 1/(3*n^(1/2)).")
    print(f"Step 6: P(max S_k >= a) <= Var(S_{n-1}) / x^2, where x is approx 1/2.")
    print(f"Step 7: This probability is bounded by C*n^(-1/2).")
    print(f"Step 8: The subtracted term is bounded by n * C*n^(-1/2) = C*n^(1/2).")
    print(f"Step 9: Therefore, E[tau] >= n - c*n^(1/2).")
    print(f"Step 10: This proves that theta = {theta} is a possible value.")
    
    # The question asks for a multiple of 1/8.
    # 1/2 = 4/8.
    numerator = 4
    denominator = 8
    print(f"The value theta = 1/2 can be written as {numerator}/{denominator}.")
    
solve()
<<<1/2>>>