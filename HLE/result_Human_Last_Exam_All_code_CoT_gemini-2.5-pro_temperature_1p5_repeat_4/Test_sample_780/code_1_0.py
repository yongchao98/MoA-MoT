import sys

def solve():
    """
    Solves the problem by calculating S(N) mod p using matrix exponentiation
    and properties of modular arithmetic.
    """
    p = 23627
    
    # Parameters from the problem
    # X is the number of ways to color a 2x1 column
    # k is the number of 'bad' colors
    X = 510**2
    k = 203
    
    # Parameters modulo p
    X_mod_p = X % p
    k_mod_p = k % p
    
    # Recurrence coefficient C = X - 1 mod p
    C = (X_mod_p - 1 + p) % p
    
    # Check if X_mod_p == k_mod_p, which is the key simplification
    # assert X_mod_p == k_mod_p
    
    # Initial values for the recurrence S(n)
    # S(1) = X
    # S(2) = X^2
    s1 = X_mod_p
    s2 = (X_mod_p**2) % p
    
    # Transition Matrix T
    # T = [[C, C], [1, 0]]
    
    # We need to compute T^(-2) * [s2, s1]^T
    # First, find T^(-1)
    
    # Determinant of T
    det_T = (C * 0 - C * 1 + p) % p
    
    # Modular inverse of det_T
    # pow(x, -1, m) is available from Python 3.8
    if sys.version_info.major == 3 and sys.version_info.minor >= 8:
        inv_det_T = pow(det_T, -1, p)
    else:
        # For older python versions, we implement the extended Euclidean algorithm
        def extended_gcd(a, b):
            if a == 0:
                return b, 0, 1
            d, x1, y1 = extended_gcd(b % a, a)
            x = y1 - (b // a) * x1
            y = x1
            return d, x, y

        d, x, y = extended_gcd(det_T, p)
        inv_det_T = x % p

    
    # Adjugate matrix of T is [[0, -C], [-1, C]]
    # T^(-1) = inv_det_T * Adjugate(T)
    t_inv_00 = (inv_det_T * 0) % p
    t_inv_01 = (inv_det_T * (-C + p)) % p
    t_inv_10 = (inv_det_T * (-1 + p)) % p
    t_inv_11 = (inv_det_T * C) % p
    
    # T^(-2) = T^(-1) * T^(-1)
    t_inv2_00 = (t_inv_00 * t_inv_00 + t_inv_01 * t_inv_10) % p
    t_inv2_01 = (t_inv_00 * t_inv_01 + t_inv_01 * t_inv_11) % p
    t_inv2_10 = (t_inv_10 * t_inv_00 + t_inv_11 * t_inv_10) % p
    t_inv2_11 = (t_inv_10 * t_inv_01 + t_inv_11 * t_inv_11) % p
    
    # Result vector [S(N), S(N-1)]^T = T^(-2) * [s2, s1]^T
    sN = (t_inv2_00 * s2 + t_inv2_01 * s1) % p
    
    print(f"The modulus is p = {p}.")
    print(f"The number of ways to color a column is X = 510^2, which is {X_mod_p} mod {p}.")
    print(f"The number of 'bad' colors is k = {k}.")
    print(f"The recurrence S(n) mod p is S(n) = {C}*S(n-1) + {C}*S(n-2).")
    print(f"Initial values: S(1) = {s1}, S(2) = {s2}.")
    print(f"The problem reduces to computing T^(-2) * [S(2), S(1)]^T mod {p}.")
    print(f"T^(-1) = [[{t_inv_00}, {t_inv_01}], [{t_inv_10}, {t_inv_11}]] mod {p}.")
    print(f"T^(-2) = [[{t_inv2_00}, {t_inv2_01}], [{t_inv2_10}, {t_inv2_11}]] mod {p}.")
    print(f"The final calculation is S(N) mod {p} = ({t_inv2_00} * {s2} + {t_inv2_01} * {s1}) mod {p}.")
    print(f"The result is {sN}.")

solve()
<<<6575>>>