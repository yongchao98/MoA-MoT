def solve():
    """
    Calculates the number of non-zero terms in the asymptotic expansion of f(x)
    up to and including the term in x^{-100}.
    The coefficients c_k of the expansion f(x) ~ sum_{k=2 to inf} c_k x^{-k}
    are determined by the following recurrence relations:
    - c_2 = 1
    - c_k = 1 for odd k >= 3
    - c_k = 1 - c_{k/2} for even k >= 4
    """
    # Using a dictionary for memoization of coefficients.
    coeffs = {}

    # Calculate coefficients c_k for k from 2 to 100.
    for k in range(2, 101):
        if k % 2 == 1:
            # For odd k, c_k is 1.
            coeffs[k] = 1
        else: # k is even
            if k == 2:
                # Base case for the recurrence.
                coeffs[k] = 1
            else:
                # For even k > 2, c_k = 1 - c_{k/2}.
                # The int() conversion is implicit in Python 2 for k/2,
                # but explicit k//2 is safer in Python 3.
                coeffs[k] = 1 - coeffs[k // 2]
    
    # Count the number of non-zero coefficients.
    nonzero_count = 0
    for k in range(2, 101):
        if coeffs[k] != 0:
            nonzero_count += 1
            
    print(nonzero_count)

solve()