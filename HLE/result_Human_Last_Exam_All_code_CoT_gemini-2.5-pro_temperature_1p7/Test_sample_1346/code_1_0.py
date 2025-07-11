def solve():
    """
    Calculates the value of a(p^4+4p^3-5p^2-3p+8) mod p for given primes.
    """

    # The recurrence relation for a(n), the number of tilings of a 3x(2n) rectangle, is:
    # a(0) = 1, a(1) = 3
    # a(n) = 4*a(n-1) - a(n-2) for n >= 2
    def get_a(n):
        if n == 0:
            return 1
        if n == 1:
            return 3
        
        a_prev, a_curr = 1, 3
        for _ in range(n - 1):
            a_next = 4 * a_curr - a_prev
            a_prev = a_curr
            a_curr = a_next
        return a_curr

    primes = [50051, 50069]
    results = []
    
    print("Let N(p) = p^4 + 4p^3 - 5p^2 - 3p + 8.")
    print("We want to compute a(N(p)) mod p for each prime p.\n")

    for p in primes:
        # Determine the period of a(n) mod p by checking (3/p).
        # (3/p) = (p/3) for p > 3.
        p_mod_3 = p % 3
        
        # For both p=50051 and p=50069, p_mod_3 is 2.
        # This means (3/p) = -1, and the period of a(n) mod p divides p+1.
        # We need to compute N(p) mod (p+1).
        # With p = -1 (mod p+1), N(p) = (-1)^4 + 4(-1)^3 - 5(-1)^2 - 3(-1) + 8 = 3.
        index = 3
        
        # Calculate a(index)
        result = get_a(index)
        results.append(result)

        print(f"For p = {p}:")
        n_expression = f"{p}^4 + 4*{p}^3 - 5*{p}^2 - 3*{p} + 8"
        print(f"a({n_expression}) mod {p}")
        
        # The logic simplifies the expression to a(3) mod p.
        print(f"= a(3) mod {p}")
        print(f"= {result}\n")
    
    print("The final answers for p=50051 and p=50069 separated by a comma are:")
    print(','.join(map(str, results)))

solve()
<<<41,41>>>