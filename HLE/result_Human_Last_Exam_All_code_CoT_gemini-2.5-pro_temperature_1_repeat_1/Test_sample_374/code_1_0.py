import math

def get_prime_factorization(num):
    """
    Computes the prime factorization of a given number.
    Returns a dictionary where keys are prime factors and values are their powers.
    """
    factors = {}
    # Count the number of factors of 2
    while num % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        num = num // 2
    # Count factors of odd primes
    d = 3
    while d * d <= num:
        while num % d == 0:
            factors[d] = factors.get(d, 0) + 1
            num = num // d
        d += 2
    if num > 1:
        factors[num] = factors.get(num, 0) + 1
    return factors

def main():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    n = 4  # Dimension of the vector space, from D being elementary abelian of order 16=2^4
    q = 2  # Characteristic of the field k

    # Step 1: Calculate the order of GL(n, q)
    # |GL(n, q)| = (q^n - 1)(q^n - q)...(q^n - q^(n-1))
    order_gl = 1
    terms = []
    for i in range(n):
        term = q**n - q**i
        terms.append(term)
        order_gl *= term
    
    print(f"The defect group D is elementary abelian of order 16, so D is isomorphic to a 4-dimensional vector space over F_2.")
    print(f"The automorphism group Aut(D) is isomorphic to GL(4, 2).")
    print(f"The order of GL(4, 2) is calculated as (2^4-1) * (2^4-2) * (2^4-4) * (2^4-8).")
    print(f"This is {terms[0]} * {terms[1]} * {terms[2]} * {terms[3]} = {order_gl}.")
    
    # Step 2: The order of E must be odd, so we find the largest odd divisor of |GL(4, 2)|
    odd_part = order_gl
    while odd_part % 2 == 0:
        odd_part = odd_part // 2
        
    print(f"\nThe inertial quotient E is a 2'-group, so its order must be odd.")
    print(f"The highest possible order for E is the largest odd divisor of {order_gl}, which is {odd_part}.")

    # Step 3: Find the prime factorization of the odd part to display the final equation
    odd_factors = get_prime_factorization(odd_part)
    
    equation_parts = []
    for p, exp in sorted(odd_factors.items()):
        if exp > 1:
            equation_parts.append(f"{p}**{exp}")
        else:
            equation_parts.append(str(p))
            
    equation_str = " * ".join(equation_parts)
    
    print("\nThe final calculation is:")
    print(f"{equation_str} = {odd_part}")
    
if __name__ == "__main__":
    main()