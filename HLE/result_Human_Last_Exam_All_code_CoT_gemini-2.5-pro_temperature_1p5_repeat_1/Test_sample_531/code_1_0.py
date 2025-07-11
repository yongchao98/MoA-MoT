def solve():
    """
    This function calculates the value of (m+M) mod 65539 based on the derived formulas.
    """
    p = 65539

    # Calculate S = 2^32 mod p
    # S = 2^32 = (2^16)^2 = 65536^2.
    # 65536 = 65539 - 3, so 65536 is congruent to -3 mod 65539.
    # S is congruent to (-3)^2 = 9 mod 65539.
    S_mod = pow(2, 32, p)

    # Calculate m mod p
    # m = 4*S^2 - 16*S
    m_mod = (4 * pow(S_mod, 2, p) - 16 * S_mod) % p

    # Calculate M mod p
    # M = S * (S-4)^2 * (S+8) / 27
    # We need the modular inverse of 27 modulo p.
    inv27 = pow(27, -1, p)
    
    # Calculate each term of M modulo p
    term1 = S_mod
    term2 = pow(S_mod - 4, 2, p)
    term3 = (S_mod + 8) % p
    
    # Combine the terms
    M_mod = (term1 * term2 * term3 * inv27) % p
    
    # Calculate the final result (m+M) mod p
    result = (m_mod + M_mod) % p
    
    print(f"The second smallest value, m, modulo {p} is: {m_mod}")
    print(f"The second largest value, M, modulo {p} is: {M_mod}")
    print(f"The sum (m + M) modulo {p} is: ({m_mod} + {M_mod}) mod {p} = {result}")

solve()