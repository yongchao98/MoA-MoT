import math

def prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
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

def solve_inertial_quotient_order():
    """
    Calculates the highest possible order for the inertial quotient E based on the problem's conditions.
    """
    print("Step 1: Laying out the theoretical foundation from Block Theory.")
    print("--------------------------------------------------------------")
    print("Let B be a block of a group algebra kG, with k being a field of characteristic p.")
    print("Let D be the defect group of B and E be its inertial quotient.")
    print("A key result states that E is isomorphic to a p'-subgroup of Out(D) = Aut(D)/Inn(D).")
    print("\nIn this specific problem:")
    print("- The characteristic is p = 2, so the order of E must be odd.")
    print("- The defect group D is elementary abelian of order 16, so D ≅ (Z/2Z)^4.")
    print("- Viewing D as a 4D vector space over F_2, its automorphism group is Aut(D) ≅ GL(4, 2).")
    print("- Since D is abelian, Inn(D) is trivial, so Out(D) ≅ GL(4, 2).")
    print("\nConclusion: The highest possible order for E is the largest odd divisor of |GL(4, 2)|.")
    print("--------------------------------------------------------------\n")

    print("Step 2: Calculating the order of GL(4, 2).")
    print("--------------------------------------------------------------")
    n = 4
    q = 2
    order_gl = 1
    
    print(f"The order of GL(n, q) is given by the formula: Π[i=0 to n-1] (q^n - q^i)")
    print(f"For GL({n}, {q}), this is: |GL({n}, {q})| = (2^4 - 2^0) * (2^4 - 2^1) * (2^4 - 2^2) * (2^4 - 2^3)")
    print(f"= ({16-1}) * ({16-2}) * ({16-4}) * ({16-8})")
    
    calc_terms = []
    for i in range(n):
        term = q**n - q**i
        order_gl *= term
        calc_terms.append(str(term))
        
    print(f"= {' * '.join(calc_terms)}")
    print(f"= {order_gl}")
    print("--------------------------------------------------------------\n")

    print("Step 3: Finding the largest odd divisor of |GL(4, 2)|.")
    print("--------------------------------------------------------------")
    factors = prime_factorization(order_gl)
    
    factor_str_parts = [f"{p}^{e}" for p, e in sorted(factors.items())]
    print(f"The prime factorization of {order_gl} is: {' * '.join(factor_str_parts)}")
    
    largest_odd_divisor = 1
    odd_factors_str = []
    values_str = []
    
    # Sort by prime base for consistent output
    for p, e in sorted(factors.items()):
        if p != 2:
            term_val = p**e
            largest_odd_divisor *= term_val
            values_str.append(str(term_val))
            if e > 1:
                odd_factors_str.append(f"{p}^{e}")
            else:
                odd_factors_str.append(f"{p}")

    equation = " * ".join(odd_factors_str)
    values = " * ".join(values_str)
    
    print("\nTo get the largest odd divisor, we take the product of all prime factors other than 2.")
    print(f"Final Calculation: {equation} = {values} = {largest_odd_divisor}")
    print("\nThus, the highest possible order for the inertial quotient E is", largest_odd_divisor)
    print("--------------------------------------------------------------")

# Execute the main function
solve_inertial_quotient_order()