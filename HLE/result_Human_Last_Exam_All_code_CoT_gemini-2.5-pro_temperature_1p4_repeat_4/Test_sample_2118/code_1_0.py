def count_nonzero_terms():
    """
    Calculates the number of nonzero terms in the asymptotic expansion of f(x)
    up to and including the term in x^-100.
    """
    # c[k] will store the coefficient of the x^{-k} term.
    # The sum for f(x) starts at k=2, so we can consider c[1] = 0.
    c = {1: 0}
    
    nonzero_odd_count = 0
    nonzero_even_count = 0

    # We need to find the number of non-zero terms for k from 2 to 100.
    for k in range(2, 101):
        if k % 2 == 1:
            # For odd k, c_k = 1.
            c[k] = 1
        else:
            # For even k, c_k = 1 - c_{k/2}.
            c[k] = 1 - c[k // 2]
        
        # Count non-zero coefficients, separating odd and even k
        if c[k] != 0:
            if k % 2 == 1:
                nonzero_odd_count += 1
            else:
                nonzero_even_count += 1
                
    total_nonzero_count = nonzero_odd_count + nonzero_even_count

    print("Number of non-zero terms with odd powers: {}".format(nonzero_odd_count))
    print("Number of non-zero terms with even powers: {}".format(nonzero_even_count))
    print("Total number of non-zero terms up to x^-100 is: {} + {} = {}".format(nonzero_odd_count, nonzero_even_count, total_nonzero_count))

count_nonzero_terms()
