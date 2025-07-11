def solve():
    """
    Calculates (m+M) mod 65539 based on the derived formulas.
    """
    p = 65539

    # S = 2^32. We only need its value modulo p.
    # S = 2^32 = (2^16)^2 = 65536^2
    # 65536 = 65539 - 3, so 65536 = -3 (mod p)
    # S = (-3)^2 = 9 (mod p)
    S_mod_p = 9

    # Calculate m mod p
    # m = 4 * S * (S - 4)
    m_mod_p = (4 * S_mod_p * (S_mod_p - 4)) % p

    # Calculate M mod p
    # M = S * (S-4)^2 * (S+8) / 27
    # We need the modular inverse of 27 mod p.
    # pow(27, -1, p) in Python 3.8+ or pow(27, p-2, p) by Fermat's Little Theorem
    inv_27 = pow(27, -1, p)
    
    # Numerator of M
    num_M_mod_p = (S_mod_p * pow(S_mod_p - 4, 2, p) * (S_mod_p + 8)) % p
    
    # M mod p
    M_mod_p = (num_M_mod_p * inv_27) % p

    # Final result
    result = (m_mod_p + M_mod_p) % p
    
    print(f"The second smallest value, m, modulo {p} is: {m_mod_p}")
    print(f"The second largest value, M, modulo {p} is: {M_mod_p}")
    print(f"The final equation is ({m_mod_p} + {M_mod_p}) mod {p}")
    print(f"The result is: {result}")

solve()
<<<22668>>>