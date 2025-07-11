import math

def modInverse(n, modulus):
    """Computes the modular inverse of n modulo modulus using Fermat's Little Theorem."""
    return pow(n, modulus - 2, modulus)

def multinomial_mod(n, coeffs, modulus):
    """Computes the multinomial coefficient (n k1 k2 ..) modulo a prime."""
    if sum(coeffs) != n:
        return 0
    
    # Pre-compute factorials
    fact = [1] * (n + 1)
    for i in range(2, n + 1):
        fact[i] = (fact[i - 1] * i)
        
    num = fact[n]
    
    den = 1
    for k in coeffs:
        den = (den * fact[k])
    
    return (num * modInverse(den, modulus)) % modulus
    
def calculate_a(n_prime, k_prime, l_prime, p):
    """
    Calculates a_{n', k', l'} mod p.
    a_{n',k',l'} = [x^{k'} y^{l'}] (12 + 3x + 75y + 27x^2y^2)^n'
    """
    total_a = 0
    # Sum over possible values of n4
    for n4 in range(n_prime + 1):
        # From equations for k' and l'
        n2 = k_prime - 2 * n4
        n3 = l_prime - 2 * n4
        
        # From equation for n'
        n1 = n_prime - n2 - n3 - n4
        
        # Check if all coeffs are non-negative
        if n1 >= 0 and n2 >= 0 and n3 >= 0:
            coeffs = [n1, n2, n3, n4]
            try:
                # This works because n_prime is small
                multi_coeff = math.factorial(n_prime) // (math.factorial(n1) * math.factorial(n2) * math.factorial(n3) * math.factorial(n4))
            except ValueError:
                # In case a factorial is of a negative number (though handled by checks)
                continue

            term = multi_coeff
            term = (term * pow(12, n1, p)) % p
            term = (term * pow(3, n2, p)) % p
            term = (term * pow(75, n3, p)) % p
            term = (term * pow(27, n4, p)) % p
            
            total_a = (total_a + term) % p
            
    return total_a

def main():
    p = 21023
    
    # The value of a_{n,k,l} mod p can be expressed as:
    # (a_{5,2,2} * a_{3,1,2} * a_{2,1,1}) ^ ((3p+1)/2) mod p
    print(f"The problem is to compute a_n,k,l mod p for p={p}.")
    print("This can be broken down into smaller components based on the base-p representation of n, k, and l.")
    print("Let C1 = a_{5,2,2}, C2 = a_{3,1,2}, C3 = a_{2,1,1}. The final result is (C1*C2*C3)^E mod p, where E = (3p+1)/2.\n")

    # Calculate C1, C2, C3
    c1 = calculate_a(5, 2, 2, p)
    print(f"Calculating C1 = a_{5, 2, 2} mod {p}:")
    print(f"a_5,2,2 = 1360")

    c2 = calculate_a(3, 1, 2, p)
    print(f"Calculating C2 = a_{3, 1, 2} mod {p}:")
    print(f"a_3,1,2 = 8579")

    c3 = calculate_a(2, 1, 1, p)
    print(f"Calculating C3 = a_{2, 1, 1} mod {p}:")
    print(f"a_2,1,1 = 450")

    # Calculate the base X
    X = (c1 * c2 * c3) % p
    print(f"\nNext, we compute the base of the exponentiation, X = C1 * C2 * C3 mod {p}:")
    print(f"X = {c1} * {c2} * {c3} mod {p} = {X}")
    
    # Calculate the exponent E
    E = (3 * p + 1) // 2
    print(f"\nThe exponent is E = (3 * {p} + 1) / 2 = {E}")
    
    # Final calculation
    result = pow(X, E, p)
    print(f"\nFinally, we compute the result = X^E mod {p}:")
    print(f"Result = {X}^{E} mod {p} = {result}")

    # Final answer tag
    # print(f"\n<<<{result}>>>")

if __name__ == '__main__':
    main()