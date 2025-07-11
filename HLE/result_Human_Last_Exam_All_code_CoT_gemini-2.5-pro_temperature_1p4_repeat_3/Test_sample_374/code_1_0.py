import math

def get_prime_factorization(n):
    """Computes the prime factorization of n as a dictionary."""
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def solve_problem():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    p = 2
    n = 4
    order_D = p**n

    print(f"The problem states that the defect group D has order {order_D} and is elementary abelian.")
    print(f"This means D is isomorphic to a {n}-dimensional vector space over the field with {p} elements, F_{p}.")
    print("The inertial quotient E is a p'-group that acts faithfully on D.")
    print(f"Therefore, E is isomorphic to an odd-order subgroup of Aut(D) = GL({n}, {p}).")
    print("The highest possible order for E is the odd part of the order of GL(4, 2).\n")

    # Calculate the order of GL(n, p)
    order_gl = 1
    terms = []
    print("The order of GL(n, p) is given by the formula: (p^n - 1) * (p^n - p) * ... * (p^n - p^(n-1))")
    print(f"For GL({n}, {p}), this is:")
    
    equation_parts = []
    for i in range(n):
        term = (p**n - p**i)
        terms.append(term)
        order_gl *= term
        equation_parts.append(str(term))

    print(f"|GL({n}, {p})| = {' * '.join(equation_parts)}")
    print(f"|GL({n}, {p})| = {order_gl}\n")

    # Find the odd part of the order
    odd_part = order_gl
    while odd_part % p == 0:
        odd_part //= p

    print(f"The order of GL({n}, {p}) is {order_gl}.")
    print(f"The characteristic is p = {p}, so we need to find the largest factor of {order_gl} that is not divisible by {p}.")
    print(f"This is the odd part of the order, which is {odd_part}.\n")
    
    print("To show the final equation for this value, we find its prime factors.")
    # We find prime factors of the odd part to display the equation
    # We know p=2, so we only need to factor odd numbers
    factors_map = {}
    temp_n_odd = odd_part
    d = 3
    while d * d <= temp_n_odd:
        while temp_n_odd % d == 0:
            factors_map[d] = factors_map.get(d, 0) + 1
            temp_n_odd //= d
        d += 2
    if temp_n_odd > 1:
        factors_map[temp_n_odd] = factors_map.get(temp_n_odd, 0) + 1
    
    number_terms = []
    for factor, power in sorted(factors_map.items()):
        number_terms.append(str(int(factor**power)))

    equation_str = " * ".join(number_terms)

    print(f"The highest possible order for E is {odd_part}.")
    print("The final equation is:")
    print(f"{equation_str} = {odd_part}")


solve_problem()

# The final numerical answer is the odd part of the order calculated.
final_answer = 315
print(f"\n<<<{final_answer}>>>")