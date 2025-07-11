def count_nonzero_terms():
    """
    Determines the number of nonzero terms in the asymptotic expansion of f(x)
    up to and including the term in x^-100.
    """
    limit = 100
    nonzero_exponents = []

    for m in range(2, limit + 1):
        # For each integer m, find its unique representation as m = q * 2^v,
        # where q is odd.
        temp_m = m
        v = 0
        while temp_m > 0 and temp_m % 2 == 0:
            v += 1
            temp_m //= 2
        q = temp_m

        is_nonzero = False
        if q == 1:
            # Case 1: m is a power of 2, m = 2^v.
            # The coefficient is non-zero if v is odd.
            if v % 2 != 0:
                is_nonzero = True
        else:
            # Case 2: m is not a power of 2, m = q * 2^v with q > 1 odd.
            # The coefficient is non-zero if v is even.
            if v % 2 == 0:
                is_nonzero = True
        
        if is_nonzero:
            nonzero_exponents.append(m)
            
    print("The exponents of the non-zero terms up to x^-100 are:")
    print(nonzero_exponents)
    print("\nThe total number of non-zero terms is:")
    print(len(nonzero_exponents))

count_nonzero_terms()