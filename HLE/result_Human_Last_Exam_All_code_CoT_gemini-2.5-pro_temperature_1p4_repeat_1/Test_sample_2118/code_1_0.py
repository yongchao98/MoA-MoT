def count_nonzero_terms(limit):
    """
    Calculates the number of nonzero terms in the asymptotic expansion up to x^-(limit).
    """

    # Case 1: k is odd.
    # a_k = 1 for all odd k >= 3.
    # We count odd numbers from 3 to limit.
    count_odd = (limit - 1) // 2  # This covers 3, 5, ..., limit or limit-1

    # Case 2: k is a power of 2, k = 2^p.
    # a_k is non-zero if p is odd.
    count_even_m1 = 0
    p = 1
    while True:
        k = 2**p
        if k > limit:
            break
        if p % 2 != 0:  # p must be odd
            count_even_m1 += 1
        p += 1

    # Case 3: k = m * 2^p where m >= 3 is odd and p >= 1.
    # a_k is non-zero if p is even.
    count_even_m_ge_3 = 0
    p = 2
    while True:
        power_of_2 = 2**p
        if power_of_2 > limit:
            break
        if p % 2 == 0: # p must be even
            m = 3
            while True:
                k = m * power_of_2
                if k > limit:
                    break
                count_even_m_ge_3 += 1
                m += 2 # Next odd m
        p += 1
        
    total_count = count_odd + count_even_m1 + count_even_m_ge_3
    
    # Print the breakdown of the calculation as requested
    print(f"Number of non-zero terms for odd k: {count_odd}")
    print(f"Number of non-zero terms for k=2^p: {count_even_m1}")
    print(f"Number of non-zero terms for other even k: {count_even_m_ge_3}")
    print(f"\nThe total number of non-zero terms is calculated as:")
    print(f"{count_odd} + {count_even_m1} + {count_even_m_ge_3} = {total_count}")
    return total_count

if __name__ == "__main__":
    limit = 100
    count_nonzero_terms(limit)
