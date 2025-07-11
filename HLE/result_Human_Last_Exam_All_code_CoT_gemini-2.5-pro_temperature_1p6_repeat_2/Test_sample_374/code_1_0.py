import math
from collections import Counter

def get_prime_factorization(n):
    """Returns a dictionary with prime factors of n and their powers."""
    factors = Counter()
    d = 2
    temp = n
    # Find factors up to sqrt(n)
    while d * d <= temp:
        if temp % d == 0:
            count = 0
            while temp % d == 0:
                count += 1
                temp //= d
            factors[d] = count
        d += 1
    # If a prime factor greater than sqrt(n) remains
    if temp > 1:
       factors[temp] = 1
    return factors

def solve():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    p = 2
    n = 4 # Since the defect group D has order 16=2^4 and is elementary abelian.

    print("The problem asks for the highest possible order of the inertial quotient E of a block B.")
    print("The given information is:")
    print(f"- The characteristic of the field k is p = {p}.")
    print(f"- The defect group D of B is elementary abelian of order 16. So, D is isomorphic to (C_2)^4.")
    print("-" * 50)
    
    print("Step 1: Relate the inertial quotient E to the defect group D.")
    print("By theory, the inertial quotient E is a p'-group that is isomorphic to a subgroup of Out(D) = Aut(D)/Inn(D).")
    print("Since D is abelian, its inner automorphism group Inn(D) is trivial. Therefore, Out(D) is equal to Aut(D).")
    print(f"For an elementary abelian group D = (C_p)^n, Aut(D) is the general linear group GL(n, p).")
    print(f"In this case, Out(D) is isomorphic to GL({n}, {p}).")
    print("-" * 50)

    print(f"Step 2: Calculate the order of GL({n}, {p}).")
    print(f"The order of GL(n, p) is given by the product: (p^n - 1) * (p^n - p) * ... * (p^n - p^(n-1))")
    
    terms = []
    order_gl = 1
    for i in range(n):
        term = (p**n - p**i)
        terms.append(term)
        order_gl *= term
        
    term_str = " * ".join(map(str, terms))
    print(f"|GL({n}, {p})| = (2^4 - 1) * (2^4 - 2) * (2^4 - 4) * (2^4 - 8)")
    print(f"|GL({n}, {p})| = {term_str} = {order_gl}")
    print("-" * 50)

    print("Step 3: Find the highest possible order for E.")
    print("E must be a p'-group. Since p=2, the order of E must be odd.")
    print("The highest possible order for E is the largest odd divisor of the order of GL(n, p).")
    
    odd_order = order_gl
    power_of_p = 0
    # Remove all factors of p from the total order
    while odd_order > 0 and odd_order % p == 0:
        odd_order //= p
        power_of_p += 1

    print(f"The full order of GL({n}, {p}) is {order_gl}.")
    print(f"Factoring out the powers of {p}: {order_gl} = {p}^{power_of_p} * {odd_order}.")
    print(f"The highest possible order for E is the odd part, which is {odd_order}.")
    print("-" * 50)

    print("Step 4: Display the final calculation.")
    factors = get_prime_factorization(odd_order)
    
    num_strings = []
    for prime, power in sorted(factors.items()):
        num_strings.append(str(prime**power))

    num_equation_lhs = " * ".join(num_strings)

    print("The final result is calculated from the prime factors of the odd part of the order:")
    print(f"{num_equation_lhs} = {odd_order}")

if __name__ == '__main__':
    solve()