def solve_asymptotic_expansion():
    """
    Determines the number of non-zero terms in the asymptotic expansion of f(x)
    up to the term in x^-100.
    """
    
    # a stores the coefficients a_k of the x^-k term.
    # We know the expansion starts at k=2, so a_1=0.
    a = {1: 0}
    
    # These lists will store the indices k for which a_k is non-zero.
    nonzero_odd_indices = []
    nonzero_even_indices = []

    # Loop through powers from k=2 to k=100.
    for k in range(2, 101):
        # For odd k, a_k = 1.
        if k % 2 == 1:
            a[k] = 1
        # For even k, a_k = 1 - a_{k/2}.
        else:
            a[k] = 1 - a[k // 2]
        
        # Check if the computed coefficient is non-zero and categorize it.
        if a[k] != 0:
            if k % 2 == 1:
                nonzero_odd_indices.append(k)
            else:
                nonzero_even_indices.append(k)

    # Count the number of terms in each category.
    count_odd = len(nonzero_odd_indices)
    count_even = len(nonzero_even_indices)
    total_count = count_odd + count_even

    # Output the result, showing the breakdown as requested.
    print(f"Number of non-zero terms from odd powers (x^-(2n+1)) up to x^-100: {count_odd}")
    print(f"Number of non-zero terms from even powers (x^-(2n)) up to x^-100: {count_even}")
    print(f"Total number of non-zero terms = {count_odd} + {count_even} = {total_count}")

solve_asymptotic_expansion()