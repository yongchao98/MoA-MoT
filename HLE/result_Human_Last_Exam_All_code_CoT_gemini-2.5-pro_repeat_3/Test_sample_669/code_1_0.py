import sys

# It might take a few seconds to run

def solve():
    """
    Solves the problem by finding the value of a_n,k,l mod p.
    """
    p = 21023

    # A dictionary to cache factorials and their modular inverses to speed up computation.
    fact_cache = {0: 1}
    inv_fact_cache = {0: 1}

    def mod_inverse(n, mod):
        """Computes modular inverse of n modulo mod using Fermat's Little Theorem."""
        return pow(n, mod - 2, mod)

    def factorial(n, mod):
        """Computes n! mod p."""
        if n in fact_cache:
            return fact_cache[n]
        res = fact_cache[max(fact_cache.keys())]
        for i in range(max(fact_cache.keys()) + 1, n + 1):
            res = (res * i) % mod
            fact_cache[i] = res
        return res

    def inv_factorial(n, mod):
        """Computes (n!)^-1 mod p."""
        if n in inv_fact_cache:
            return inv_fact_cache[n]
        # Ensure factorial(n, mod) is computed and cached first
        fact_n = factorial(n, mod)
        inv_fact_n = mod_inverse(fact_n, mod)
        inv_fact_cache[n] = inv_fact_n
        return inv_fact_n

    def multinomial_coeff(n, ks, mod):
        """Computes the multinomial coefficient n! / (k1! * k2! * ...) mod p."""
        if sum(ks) != n:
            return 0
        
        numerator = factorial(n, mod)
        denominator = 1
        for k in ks:
            if k < 0: return 0
            denominator = (denominator * factorial(k, mod)) % mod
            
        return (numerator * mod_inverse(denominator, mod)) % mod

    def c(n, k, l, mod):
        """
        Calculates [x^k y^l] (12 + 3x + 75y + 27x^2y^2)^n mod p.
        This is based on the multinomial expansion:
        Sum over n0+n1+n2+n3=n of (n!/(n0!n1!n2!n3!)) * 12^n0 * 3^n1 * 75^n2 * 27^n3
        where the powers of x and y must match k and l:
        n1 + 2*n3 = k
        n2 + 2*n3 = l
        """
        total = 0
        # Iterate over possible values of n3
        for n3 in range(min(k // 2, l // 2) + 1):
            # From power constraints
            n1 = k - 2 * n3
            n2 = l - 2 * n3
            
            # From sum constraint
            n0 = n - n1 - n2 - n3
            
            if n0 < 0 or n1 < 0 or n2 < 0:
                continue
            
            # Calculate multinomial coefficient
            mult_coeff = multinomial_coeff(n, [n0, n1, n2, n3], mod)
            
            # Calculate powers of coefficients
            term_val = pow(12, n0, mod)
            term_val = (term_val * pow(3, n1, mod)) % mod
            term_val = (term_val * pow(75, n2, mod)) % mod
            term_val = (term_val * pow(27, n3, mod)) % mod
            
            term = (mult_coeff * term_val) % mod
            total = (total + term) % mod
            
        return total

    # The base-p digits of (n,k,l) repeat in a 3-cycle:
    # d0=(5,2,2), d1=(3,1,2), d2=(2,1,1)
    
    print(f"The prime modulus is p = {p}")
    print("Calculating coefficients for one block of digits...")

    # Calculate coefficient for the first set of digits (n_j, k_j, l_j) = (5, 2, 2)
    C0 = c(5, 2, 2, p)
    print(f"c(n=5, k=2, l=2) mod {p} = {C0}")
    
    # Calculate coefficient for the second set of digits (n_j, k_j, l_j) = (3, 1, 2)
    C1 = c(3, 1, 2, p)
    print(f"c(n=3, k=1, l=2) mod {p} = {C1}")

    # Calculate coefficient for the third set of digits (n_j, k_j, l_j) = (2, 1, 1)
    C2 = c(2, 1, 1, p)
    print(f"c(n=2, k=1, l=1) mod {p} = {C2}")

    # The result is the product of these coefficients raised to the number of blocks
    X = (C0 * C1 * C2) % p
    print(f"The product of coefficients for one block is X = ({C0} * {C1} * {C2}) mod {p} = {X}")

    # The number of blocks is (3p+1)/2
    E = (3 * p + 1) // 2
    print(f"The exponent, which is the number of blocks, is E = (3*{p} + 1)/2 = {E}")
    
    # Final calculation using modular exponentiation
    final_result = pow(X, E, p)
    print(f"The final result is {X}^{E} mod {p} = {final_result}")
    
    return final_result

result = solve()
print(f"\nFinal answer: {result}")
<<<6881>>>