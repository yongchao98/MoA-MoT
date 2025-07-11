import math

def solve():
    """
    Calculates the value of a_{n,k,l} mod p based on the given recurrence and parameters.
    """
    p = 21023

    # Pre-compute factorials modulo p up to n=5, which is the max needed.
    factorials = [1] * 6
    for i in range(2, 6):
        factorials[i] = (factorials[i - 1] * i) % p

    def mod_inverse(n, modulus):
        """Computes modular inverse using Fermat's Little Theorem."""
        return pow(n, modulus - 2, modulus)

    def multinomial_coeff(n, k_list, modulus):
        """Computes multinomial coefficient (n choose k1, k2, ...) mod p."""
        if sum(k_list) != n:
            return 0
        
        num = factorials[n]
        den = 1
        for k in k_list:
            if k < 0: return 0 # Should not happen with current logic, but good practice
            den = (den * factorials[k]) % modulus
            
        return (num * mod_inverse(den, modulus)) % modulus

    def calculate_a(n, k, l, modulus):
        """
        Calculates a_{n,k,l} for small n, k, l by finding the coefficient of
        x^k * y^l in (12 + 3x + 75y + 27x^2y^2)^n mod p.
        """
        total_coeff = 0
        
        # Iterate over possible values of n4 in the multinomial expansion
        # (12)^n1 * (3x)^n2 * (75y)^n3 * (27x^2y^2)^n4
        # Exponents must satisfy:
        # n1 + n2 + n3 + n4 = n
        # n2 + 2*n4 = k
        # n3 + 2*n4 = l
        for n4 in range(min(k // 2, l // 2) + 1):
            n2 = k - 2 * n4
            n3 = l - 2 * n4
            n1 = n - k - l + 3 * n4
            
            if n1 >= 0:
                coeffs = [n1, n2, n3, n4]
                multi_coeff = multinomial_coeff(n, coeffs, modulus)
                
                term_val = multi_coeff
                term_val = (term_val * pow(12, n1, modulus)) % modulus
                term_val = (term_val * pow(3, n2, modulus)) % modulus
                term_val = (term_val * pow(75, n3, modulus)) % modulus
                term_val = (term_val * pow(27, n4, modulus)) % modulus
                
                total_coeff = (total_coeff + term_val) % modulus
                
        return total_coeff

    # The base-p digits of (n,k,l) are periodic with period 3.
    # We calculate the corresponding a_{n_j, k_j, l_j} for each part of the cycle.
    
    # For j mod 3 = 0: (n_j, k_j, l_j) = (5, 2, 2)
    C0 = calculate_a(5, 2, 2, p)
    
    # For j mod 3 = 1: (n_j, k_j, l_j) = (3, 1, 2)
    C1 = calculate_a(3, 1, 2, p)
    
    # For j mod 3 = 2: (n_j, k_j, l_j) = (2, 1, 1)
    C2 = calculate_a(2, 1, 1, p)
    
    # The number of repeating cycles is (3p-1)/2 + 1 = (3*21023+1)/2 = 31535
    exponent = 31535
    
    # The final result is (C0 * C1 * C2)^exponent mod p
    V = (C0 * C1 * C2) % p
    
    # The exponent for the modular exponentiation can be reduced modulo (p-1)
    # by Fermat's Little Theorem.
    mod_exp = exponent % (p - 1)
    
    final_result = pow(V, mod_exp, p)
    
    # Output the components of the final calculation as requested
    print(f"The prime modulus is p = {p}.")
    print("The base coefficients for the periodic digits are:")
    print(f"  C_0 = a(5,2,2) mod p = {C0}")
    print(f"  C_1 = a(3,1,2) mod p = {C1}")
    print(f"  C_2 = a(2,1,1) mod p = {C2}")
    print(f"Their product is V = (C_0 * C_1 * C_2) mod p = {V}.")
    print(f"The number of repeating cycles is E = {exponent}.")
    print("\nThe final equation to solve is:")
    print(f"  Result = ({C0} * {C1} * {C2})^{exponent} mod {p}")
    print(f"         = {V}^{exponent} mod {p}")
    print(f"This is computed as {V}^{{{exponent} mod {p-1}}} mod {p} = {V}^{mod_exp} mod {p}.")
    print(f"\nThe final calculated value of a_{{n,k,l}} mod p is: {final_result}")

solve()
<<<11867>>>