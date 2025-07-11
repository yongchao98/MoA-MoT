def count_nonzero_terms():
    """
    This function calculates the number of nonzero terms in the asymptotic expansion
    of f(x) up to the term in x^-100.
    
    A term c_k is nonzero if k = 2^p * q (with q odd) satisfies one of:
    1. q = 1 and p is odd.
    2. q > 1 and p is even.
    """
    max_k = 100
    
    # Case 1: k = 2^p, p is odd
    counts = []
    p = 1
    k = 2**p
    count_case1 = 0
    while k <= max_k:
        count_case1 += 1
        p += 2
        k = 2**p
    counts.append(count_case1)
    
    # Case 2: k = 2^p * q, q > 1 is odd, p is even
    p = 0
    while True:
        power_of_2 = 2**p
        if power_of_2 > max_k:
            break
        
        # We need to count odd integers q in [3, max_k / power_of_2]
        min_q = 3
        max_q_val = max_k // power_of_2
        
        count_p = 0
        if max_q_val >= min_q:
            # Adjust max_q_val to the largest odd number in the range
            if max_q_val % 2 == 0:
                effective_max_q = max_q_val - 1
            else:
                effective_max_q = max_q_val
            
            if effective_max_q >= min_q:
                count_p = (effective_max_q - min_q) // 2 + 1
        
        counts.append(count_p)
        
        p += 2

    # Print the equation for the total sum
    equation_parts = [str(c) for c in counts if c > 0]
    total = sum(counts)
    print(" + ".join(equation_parts) + f" = {total}")

count_nonzero_terms()
