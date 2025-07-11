import sys

# In Python 3, pow(base, exp, mod) can handle negative exponents for modular inverse.
# In Python 2, we would need to implement our own modular inverse function.
# This code assumes Python 3.6+ for simplicity.
if sys.version_info < (3, 8):
    def inv(base, mod):
        """Modular inverse function for Python versions older than 3.8"""
        g, x, y = egcd(base, mod)
        if g != 1:
            raise Exception('modular inverse does not exist')
        return x % mod

    def egcd(a, b):
        """Extended Euclidean Algorithm"""
        if a == 0:
            return (b, 0, 1)
        else:
            g, y, x = egcd(b % a, a)
            return (g, x - (b // a) * y, y)
    
    power = lambda base, exp, mod: pow(base, exp, mod) if exp >= 0 else inv(base, mod)

else:
    power = pow


def solve():
    """
    This function calculates the value of (m+M) mod 65539 based on the derived formulas.
    """
    P = 65539
    
    # Let S = 2^32. We first compute S mod P.
    # P = 2^16 + 3, so 2^16 = -3 (mod P)
    # S = 2^32 = (2^16)^2 = (-3)^2 = 9 (mod P)
    S_mod_P = 9

    # Calculate m = 4*S*(S-4) mod P
    m_val = (4 * S_mod_P * (S_mod_P - 4)) % P
    
    # Calculate M = S*(S+8)*(S-4)^2 / 27 mod P
    # We need to compute the modular inverse of 27.
    try:
        inv27 = power(27, -1, P)
    except Exception as e:
        print(f"Error calculating modular inverse: {e}")
        return

    # M = (S * (S+8) * (S-4)^2 * 27^-1) mod P
    term1 = S_mod_P
    term2 = (S_mod_P + 8) % P
    term3 = power(S_mod_P - 4, 2, P)
    
    M_val_num = (term1 * term2 * term3) % P
    M_val = (M_val_num * inv27) % P

    # Final result
    result = (m_val + M_val) % P
    
    print(f"The second smallest value m is {m_val} (mod {P}).")
    print(f"The second largest value M is {M_val} (mod {P}).")
    print(f"The final result (m+M) mod {P} is:")
    print(f"({m_val} + {M_val}) mod {P} = {result}")

solve()
<<<12168>>>