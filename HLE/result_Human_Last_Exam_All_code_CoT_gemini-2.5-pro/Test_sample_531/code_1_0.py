def solve():
    """
    Calculates the value of (m+M) mod 65539 based on the problem description.
    """
    P = 65539

    # Step 1: Calculate S mod P
    # S = 2**32. P = 2**16 + 3.
    # S = (2**16)**2 = (P - 3)**2 = P**2 - 6P + 9
    # S mod P is 9.
    S_mod = pow(2, 32, P)

    # Step 2: Calculate m mod P
    # m = S * 4 * (S - 4)
    m_mod = (S_mod * 4 * (S_mod - 4)) % P

    # Step 3: Calculate M mod P
    # M = S * (S-4)^2 * (S+8) / 27
    # We need the modular inverse of 27 mod P.
    # Since P is prime, we can use Fermat's Little Theorem: a^(P-2) is the inverse of a mod P.
    inv27 = pow(27, P - 2, P)
    
    s_minus_4_sq = pow(S_mod - 4, 2, P)
    s_plus_8 = S_mod + 8
    
    numerator_mod = (S_mod * s_minus_4_sq * s_plus_8) % P
    M_mod = (numerator_mod * inv27) % P

    # Step 4: Calculate the final result
    result = (m_mod + M_mod) % P

    # Print the equation with the calculated numbers
    print(f"Let P = {P}")
    print(f"m = S * 4 * (S-4)")
    print(f"m mod P = ({S_mod} * 4 * ({S_mod}-4)) mod {P} = {m_mod}")
    print(f"M = S * (S-4)^2 * (S+8) / 27")
    print(f"M mod P = ({S_mod} * ({S_mod}-4)^2 * ({S_mod}+8) * 27^(-1)) mod {P} = {M_mod}")
    print(f"\nFinal Calculation:")
    print(f"(m + M) mod {P} = ({m_mod} + {M_mod}) mod {P} = {result}")

solve()