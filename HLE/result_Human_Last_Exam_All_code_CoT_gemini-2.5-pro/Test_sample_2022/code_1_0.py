def solve():
    """
    Calculates the value of F(n) for the given primes.
    """
    primes = [80039, 80077]
    results = []
    
    for p in primes:
        # Define N based on p
        # N = p**5 + 2*p**4 - 19*p**3 - 3*p**2 + 16*p + 6
        # We don't need to compute the full value of N, just its value modulo p-1 or 2(p+1).

        # Check if 5 is a quadratic residue modulo p
        leg_symbol = pow(5, (p - 1) // 2, p)
        
        if leg_symbol == 1:
            # Case 1: p = 80039
            # The period of the sequence F(n) divides p-1.
            # N_mod_p_minus_1 = (1 + 2 - 19 - 3 + 16 + 6) % (p - 1) which is 3
            # F(N) = F(3) = S(3) = 7/8
            a = 7
            b = 8
            
            print(f"For p={p}, the roots of the characteristic polynomial are in Z_p.")
            print(f"The calculation is ({a} * pow({b}, -1, {p})) % {p}")
            
            res = (a * pow(b, -1, p)) % p
            results.append(res)

        else: # leg_symbol == p - 1
            # Case 2: p = 80077
            # We use the formula S(n) = (L_n + 2*F_n - 1) / 2^n
            # The period of L_n and F_n divides 2*(p+1).
            # We found N = m*(p+1)+7 where m is even. So N is equivalent to 7.
            # L_7 = 29, F_7 = 13
            # Numerator = L_7 + 2*F_7 - 1 = 29 + 2*13 - 1 = 54
            numerator = 54
            
            # For the denominator 2^N, we need N mod (p-1)
            # N_mod_p_minus_1 = (m_mod * 2 + 7) mod (p-1)
            # m = p^4+p^3-20p^2+17p-1
            # m_mod_p_minus_1 = (1+1-20+17-1) % (p-1) = -2
            # N_mod_p_minus_1 = (-2 * 2 + 7) % (p-1) = 3
            # Denominator = 2^N mod p = 2^3 = 8
            denominator = 8
            
            print(f"For p={p}, the roots are in the extension field Z_p^2.")
            print(f"The calculation is ({numerator} * pow({denominator}, -1, {p})) % {p}")

            res = (numerator * pow(denominator, -1, p)) % p
            results.append(res)
            
    print(",".join(map(str, results)))

solve()