def count_nonzero_terms():
    """
    Calculates the number of non-zero terms in the asymptotic expansion of f(x)
    up to and including the term in x^-100.
    """
    
    # These counters correspond to the cases derived in the analysis.
    # n = 2^p * m, where m is odd.
    # a_n is non-zero if (m=1 and p is odd) or (m>=3 and p is even).
    
    # Case: m = 1, n = 2^p, p is odd
    count_m1_p_odd = 0
    
    # Case: m >= 3, p is even
    # We split this by the value of p for clarity in the output.
    count_m3_p0 = 0  # p=0, n=m
    count_m3_p2 = 0  # p=2, n=4m
    count_m3_p4 = 0  # p=4, n=16m
    count_m3_p6 = 0  # p=6, n=64m (and so on, but p=6 is the last relevant one)

    for n in range(2, 101):
        m = n
        p = 0
        while m > 0 and m % 2 == 0:
            m //= 2
            p += 1
        
        if m == 1:  # n is a power of 2
            if p % 2 != 0:
                count_m1_p_odd += 1
        else:  # m is odd and >= 3
            if p % 2 == 0:
                if p == 0:
                    count_m3_p0 += 1
                elif p == 2:
                    count_m3_p2 += 1
                elif p == 4:
                    count_m3_p4 += 1
                elif p == 6:
                    count_m3_p6 += 1
    
    total_count = count_m1_p_odd + count_m3_p0 + count_m3_p2 + count_m3_p4 + count_m3_p6
    
    print(f"The number of non-zero terms is the sum of counts from several cases based on n = 2^p * m:")
    print(f"1. n = 2^p (m=1), p is odd: {count_m1_p_odd} terms")
    print(f"2. n = m (p=0), m is odd >= 3: {count_m3_p0} terms")
    print(f"3. n = 4m (p=2), m is odd >= 3: {count_m3_p2} terms")
    print(f"4. n = 16m (p=4), m is odd >= 3: {count_m3_p4} terms")
    print(f"5. n = 64m (p=6), m is odd >= 3: {count_m3_p6} terms")
    print(f"\nThe final equation for the total count is:")
    print(f"{count_m1_p_odd} + {count_m3_p0} + {count_m3_p2} + {count_m3_p4} + {count_m3_p6} = {total_count}")

count_nonzero_terms()