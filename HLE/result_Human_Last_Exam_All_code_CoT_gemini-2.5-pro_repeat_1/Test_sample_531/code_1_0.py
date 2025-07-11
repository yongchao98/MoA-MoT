def solve():
    """
    Calculates the value of (m+M) mod 65539 based on the derived formulas for m and M.
    """
    P = 65539
    
    # Let S = 2**32. We calculate S mod P.
    # P = 65539 = 2**16 + 3
    # So, 2**16 = -3 (mod P)
    # S = 2**32 = (2**16)**2 = (-3)**2 = 9 (mod P)
    S_mod_P = 9

    # From our analysis, m = 4*S*(S-4)
    # We calculate m mod P
    m_mod_P = (4 * S_mod_P * (S_mod_P - 4)) % P

    # From our analysis, M = S * (S+8) * (S-4)**2 / 27
    # We calculate M mod P. This requires the modular inverse of 27.
    # Using pow(27, -1, P) for python 3.8+ or Extended Euclidean Algorithm
    inv27 = pow(27, -1, P)
    
    # M = S * (S+8) * (S-4)^2 * 27^(-1) (mod P)
    term1_M = S_mod_P
    term2_M = (S_mod_P + 8) % P
    term3_M = pow(S_mod_P - 4, 2, P)
    
    M_mod_P = (term1_M * term2_M * term3_M * inv27) % P

    # The final result is (m + M) mod P
    result = (m_mod_P + M_mod_P) % P
    
    # Output the numbers in the final equation as requested
    print(f"({m_mod_P} + {M_mod_P}) % {P} = {result}")

solve()