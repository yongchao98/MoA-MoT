def count_nonzero_terms():
    """
    Counts the number of nonzero terms in the asymptotic expansion of f(x)
    for powers of x from x^-2 to x^-100.
    """
    limit = 100
    
    # Case 1: n = 2^p, with p odd.
    # We count n <= 100.
    # p=1: n=2
    # p=3: n=8
    # p=5: n=32
    # p=7: n=128 > 100
    count1 = 3
    print(f"Number of nonzero terms for n = 2^p with p odd: {count1}")

    # Case 2: n = 2^p * q, with q > 1 odd and p even.
    # We count n <= 100.
    
    # p=0: n = q, where q is an odd integer >= 3.
    # q can be 3, 5, 7, ..., 99.
    # Number of terms = (99 - 3) / 2 + 1 = 49.
    count2_p0 = (99 - 3) // 2 + 1
    print(f"Number of nonzero terms for n = q (p=0) with q > 1 odd: {count2_p0}")
    
    # p=2: n = 4 * q, where q is an odd integer >= 3.
    # 4q <= 100 => q <= 25.
    # q can be 3, 5, 7, ..., 25.
    # Number of terms = (25 - 3) / 2 + 1 = 12.
    count2_p2 = (25 - 3) // 2 + 1
    print(f"Number of nonzero terms for n = 4q (p=2) with q > 1 odd: {count2_p2}")

    # p=4: n = 16 * q, where q is an odd integer >= 3.
    # 16q <= 100 => q <= 6.25.
    # q can be 3, 5.
    # Number of terms = 2.
    count2_p4 = 2
    print(f"Number of nonzero terms for n = 16q (p=4) with q > 1 odd: {count2_p4}")

    # p=6: n = 64 * q, where q is an odd integer >= 3.
    # 64q <= 100 => q <= 1.5625. No possible values for q.
    # Higher even powers of p will also yield no results.
    
    total_count = count1 + count2_p0 + count2_p2 + count2_p4
    
    print(f"\nThe total number of nonzero terms is the sum of these counts:")
    print(f"{count1} + {count2_p0} + {count2_p2} + {count2_p4} = {total_count}")
    
count_nonzero_terms()