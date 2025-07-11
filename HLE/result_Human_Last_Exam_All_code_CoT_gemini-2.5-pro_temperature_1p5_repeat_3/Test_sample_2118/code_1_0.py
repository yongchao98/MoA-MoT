def count_nonzero_terms():
    """
    This function calculates the number of non-zero terms in the asymptotic expansion
    of f(x) up to the term x^{-100}, based on the derived conditions.
    A term c_n is non-zero if:
    1. n = 2^p and p is odd.
    2. n = m * 2^p where m >= 3 is odd, and p is even.
    """
    limit = 100
    
    count_case1 = 0  # n = 2^p, p odd
    count_case2_p0 = 0 # n = m * 2^0, m >= 3 odd
    count_case2_p2 = 0 # n = m * 2^2, m >= 3 odd
    count_case2_p4 = 0 # n = m * 2^4, m >= 3 odd
    
    nonzero_exponents = []

    for n in range(2, limit + 1):
        # Decompose n into m * 2^p
        temp_n = n
        p = 0
        while temp_n % 2 == 0 and temp_n > 1:
            temp_n //= 2
            p += 1
        m = temp_n

        is_nonzero = False
        if m == 1:  # Case 1: n is a power of 2
            if p % 2 == 1:
                is_nonzero = True
                count_case1 += 1
        elif m >= 3:  # Case 2: n has an odd factor >= 3
            if p % 2 == 0:
                is_nonzero = True
                if p == 0:
                    count_case2_p0 += 1
                elif p == 2:
                    count_case2_p2 += 1
                elif p == 4:
                    count_case2_p4 += 1
        
        if is_nonzero:
            nonzero_exponents.append(n)

    total_count = count_case1 + count_case2_p0 + count_case2_p2 + count_case2_p4

    print("The non-zero terms are for n values:")
    print(nonzero_exponents)
    print("\nBreakdown of the count:")
    print(f"Number of terms for n = 2^p, p is odd: {count_case1}")
    print(f"Number of terms for n = m*2^0 (i.e., odd n >= 3): {count_case2_p0}")
    print(f"Number of terms for n = m*2^2 (e.g., 12, 20,...): {count_case2_p2}")
    print(f"Number of terms for n = m*2^4 (e.g., 48, 80,...): {count_case2_p4}")
    
    print(f"\nThe final equation for the total count is:")
    print(f"{count_case1} + {count_case2_p0} + {count_case2_p2} + {count_case2_p4} = {total_count}")

    # The final print is just the number itself as requested by user prompt
    print(f"\nTotal number of nonzero terms: {total_count}")

count_nonzero_terms()