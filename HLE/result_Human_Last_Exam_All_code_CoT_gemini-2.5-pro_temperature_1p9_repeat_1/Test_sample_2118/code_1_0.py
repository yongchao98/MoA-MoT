def solve_asymptotic_expansion():
    """
    Calculates the number of non-zero terms in the asymptotic expansion of f(x)
    up to the x^{-100} term.
    """
    # The coefficients a_k are stored in a dictionary.
    a = {}

    # The recurrence relations for the coefficients a_k are:
    # 1. a_1 = 0
    # 2. For k > 1 odd, a_k = 1
    # 3. For k = 2m (even), a_2m = 1 - a_m

    # We can compute a_k for k from 1 to 100 iteratively.
    # Base case
    a[1] = 0

    # Iteratively compute a_k for k from 2 to 100
    for k in range(2, 101):
        if k % 2 == 1:
            # Case for k odd
            a[k] = 1
        else:
            # Case for k even
            # Since we iterate k from 2 upwards, a[k // 2] will always be
            # computed before a[k].
            a[k] = 1 - a[k // 2]

    # Count the number of non-zero coefficients
    nonzero_count = 0
    for k in range(1, 101):
        if a[k] != 0:
            nonzero_count += 1
            
    total_terms = 100
    zero_count = total_terms - nonzero_count

    print(f"Determining the number of non-zero terms up to and including x^-100.")
    print(f"Total number of terms considered: {total_terms}")
    print(f"Number of zero-value terms found: {zero_count}")
    print(f"Number of non-zero terms found: {nonzero_count}")
    print(f"Final calculation: {total_terms} - {zero_count} = {nonzero_count}")

solve_asymptotic_expansion()
