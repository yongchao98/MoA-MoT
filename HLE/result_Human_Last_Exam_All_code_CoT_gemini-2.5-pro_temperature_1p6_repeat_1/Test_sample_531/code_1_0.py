import sys

def solve():
    """
    Solves the problem by calculating the second smallest (m) and second largest (M)
    values of the given expression, and then finds (m+M) mod 65539.
    """
    P = 65539
    
    # Let S = 2^32. We compute S mod P.
    # S = 2^32 = (2^16)^2 = 65536^2
    # P = 65539 = 65536 + 3, so 65536 = P - 3.
    # S = (P-3)^2 = P^2 - 6P + 9.
    # S mod P = 9.
    # We can also use pow for direct computation.
    S_mod_P = pow(2, 32, P)

    # Calculate m mod P
    # m = 4 * S * (S - 4)
    # The numbers in the equation for m are 4, S, and S-4
    m_term1 = 4
    m_term2 = S_mod_P
    m_term3 = (S_mod_P - 4 + P) % P
    
    m_mod_P = (m_term1 * m_term2 * m_term3) % P
    
    print("The second smallest value m is based on the equation 4 * S * (S - 4).")
    print(f"m mod {P} = ({m_term1} * {m_term2} * {m_term3}) mod {P} = {m_mod_P}")

    # Calculate M mod P
    # M = S * (S-4)^2 * (S+8) / 27
    inv_27 = pow(27, -1, P)
    
    # The numbers in the equation for M are S, (S-4)^2, (S+8), and 1/27
    M_term1 = S_mod_P
    M_term2 = pow(S_mod_P - 4, 2, P)
    M_term3 = (S_mod_P + 8) % P
    M_term4 = inv_27
    
    M_mod_P = (M_term1 * M_term2 * M_term3 * M_term4) % P
    
    print("\nThe second largest value M is based on the equation S * (S-4)^2 * (S+8) / 27.")
    print(f"M mod {P} = ({M_term1} * {M_term2} * {M_term3} * {M_term4}) mod {P} = {M_mod_P}")

    # Calculate the final result (m+M) mod P
    result = (m_mod_P + M_mod_P) % P
    
    print(f"\nThe value of (m + M) mod {P} is calculated as:")
    print(f"({m_mod_P} + {M_mod_P}) mod {P} = {result}")

solve()