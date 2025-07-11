import sys

def solve():
    """
    Solves the problem of calculating a_{n,k,l} mod p based on the provided plan.
    """
    
    # Define the prime modulus
    p = 21023
    
    # --- Helper functions for modular arithmetic ---
    
    def power(base, exp, mod):
        """
        Computes (base^exp) % mod using modular exponentiation.
        """
        res = 1
        base %= mod
        while exp > 0:
            if exp % 2 == 1:
                res = (res * base) % mod
            base = (base * base) % mod
            exp //= 2
        return res

    def inv(n, mod):
        """
        Computes the modular multiplicative inverse of n modulo mod.
        """
        return power(n, mod - 2, mod)

    # Use a cache for factorials as they are used multiple times for small numbers
    fact_cache = {0: 1}
    def fact(n):
        if n in fact_cache:
            return fact_cache[n]
        res = (fact(n-1) * n)
        fact_cache[n] = res
        return res
        
    def multinomial_mod(n, coeffs, mod):
        """
        Computes the multinomial coefficient (n choose k1, k2, ...) mod p.
        """
        denom = 1
        for c in coeffs:
            denom = (denom * fact(c)) % mod
        return (fact(n) * inv(denom, mod)) % mod
        
    def calculate_a_mod_p(n, k, l, mod):
        """
        Calculates a_{n,k,l} mod p for small n, k, l using the multinomial expansion
        of P(x,y)^n = (12 + 3x + 75y + 27x^2*y^2)^n.
        a_{n,k,l} = sum_{n4} C(n;n1,n2,n3,n4) * 12^n1 * 3^n2 * 75^n3 * 27^n4
        where n1=n-k-l+3*n4, n2=k-2*n4, n3=l-2*n4.
        """
        total = 0
        for n4 in range(n + 1):
            n2 = k - 2 * n4
            n3 = l - 2 * n4
            if n2 < 0 or n3 < 0:
                continue
                
            n1 = n - n2 - n3 - n4
            if n1 < 0:
                continue

            if n1 + n2 + n3 + n4 != n:
                continue
            
            multi_coeff = multinomial_mod(n, [n1, n2, n3, n4], mod)

            term = multi_coeff
            term = (term * power(12, n1, mod)) % mod
            term = (term * power(3, n2, mod)) % mod
            term = (term * power(75, n3, mod)) % mod
            term = (term * power(27, n4, mod)) % mod
            
            total = (total + term) % mod
            
        return total

    # --- Main calculation steps ---
    
    # Pre-populate the factorial cache up to the maximum needed value (5)
    for i in range(6):
        fact(i)

    # Step 1: Calculate a_{n,k,l} for the three digit tuples
    v1 = calculate_a_mod_p(5, 2, 2, p)
    v2 = calculate_a_mod_p(3, 1, 2, p)
    v3 = calculate_a_mod_p(2, 1, 1, p)
    
    # Step 2: Combine the results
    X = (v1 * v2 * v3) % p
    
    # Step 3: Determine the exponent
    # Each digit pattern repeats E = (3p+1)/2 times.
    exponent = (3 * p + 1) // 2
    
    # Step 4: Calculate the final result
    result = power(X, exponent, p)
    
    # --- Output the results step-by-step ---
    
    print(f"The prime modulus is p = {p}.")
    print("The problem simplifies to calculating smaller terms based on the base-p digits of n, k, and l.")
    print("There are three repeating tuples of digits (n_j, k_j, l_j). We calculate a_n_k_l for each.")
    print(f"1. For (n,k,l)=(5,2,2): a_{5,2,2} mod {p} = {v1}")
    print(f"2. For (n,k,l)=(3,1,2): a_{3,1,2} mod {p} = {v2}")
    print(f"3. For (n,k,l)=(2,1,1): a_{2,1,1} mod {p} = {v3}")
    
    print("\nWe combine these values into a single base X:")
    print(f"X = ({v1} * {v2} * {v3}) mod {p} = {X}")
    
    print("\nEach digit pattern repeats E times, where E = (3*p + 1)/2:")
    print(f"E = (3*{p} + 1)/2 = {exponent}")
    
    print("\nThe final answer is X^E mod p:")
    print(f"Final Result = {X}^{exponent} mod {p} = {result}")

solve()
<<<11253>>>