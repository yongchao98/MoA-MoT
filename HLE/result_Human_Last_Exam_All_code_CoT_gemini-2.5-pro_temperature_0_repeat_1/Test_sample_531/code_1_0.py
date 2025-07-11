def solve():
    """
    Calculates (m+M) mod 65539 based on the derived formulas for m and M.
    """
    p = 65539

    # S = 2^32. We calculate S mod p.
    # S = 2^32 = (2^16)^2.
    # Since p = 65539 = 2^16 + 3, we have 2^16 = p - 3, so 2^16 = -3 (mod p).
    # S = (-3)^2 = 9 (mod p).
    S_mod_p = pow(2, 32, p)

    # m = 4 * S * (S - 4)
    # We calculate m mod p.
    m_mod_p = (4 * S_mod_p * (S_mod_p - 4)) % p

    # M = S * (S - 4)^2 * (S + 8) / 27
    # We calculate M mod p. This requires the modular inverse of 27.
    inv_27 = pow(27, -1, p)
    
    term1 = S_mod_p
    term2 = pow(S_mod_p - 4, 2, p)
    term3 = (S_mod_p + 8)
    
    M_mod_p = (term1 * term2 * term3 * inv_27) % p

    # Final result is (m + M) mod p
    result = (m_mod_p + M_mod_p) % p

    # Output the components of the final equation as requested.
    print(f"The modulus is p = {p}")
    print(f"The value of S = 2^32 modulo p is: {S_mod_p}")
    print(f"The second smallest value m modulo p is: {m_mod_p}")
    print(f"The second largest value M modulo p is: {M_mod_p}")
    print(f"The final equation is ({m_mod_p} + {M_mod_p}) mod {p}")
    print(f"The result is: {result}")

solve()