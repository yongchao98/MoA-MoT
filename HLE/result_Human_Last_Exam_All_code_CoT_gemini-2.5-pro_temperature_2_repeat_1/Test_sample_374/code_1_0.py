import collections

def solve():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    n = 4
    q = 2

    print("Step 1: Understand the theoretical background.")
    print("The defect group D is elementary abelian of order 16, so D is a 4-dimensional vector space over F_2.")
    print("The automorphism group Aut(D) is isomorphic to GL(4, F_2).")
    print("The inertial quotient E is a subgroup of Aut(D) of odd order (since char(k)=2).")
    print("The maximum possible order of E is the odd part of the order of GL(4, F_2).\n")

    print("Step 2: Calculate the order of GL(4, F_2).")
    print(f"The formula for the order of GL(n, q) is the product of (q^n - q^i) for i from 0 to n-1.")
    
    terms = []
    for i in range(n):
        term = q**n - q**i
        terms.append(term)
        
    equation_str = " * ".join(map(str, terms))
    total_order = 1
    for term in terms:
        total_order *= term
    print(f"|GL(4, 2)| = (2^4 - 2^0) * (2^4 - 2^1) * (2^4 - 2^2) * (2^4 - 2^3)")
    print(f"           = {equation_str} = {total_order}\n")

    print("Step 3: Find the prime factorization of the odd part of the order.")
    prime_factors = collections.defaultdict(int)
    
    for term in terms:
        odd_part = term
        while odd_part > 0 and odd_part % 2 == 0:
            odd_part //= 2
        
        # Prime factorize the odd part
        d = 3
        temp_num = odd_part
        while d * d <= temp_num:
            while temp_num % d == 0:
                prime_factors[d] += 1
                temp_num //= d
            d += 2
        if temp_num > 1:
            prime_factors[temp_num] += 1
            
    # Sort factors for consistent output
    sorted_factors = sorted(prime_factors.items())

    print("The odd part of |GL(4, 2)| is the product of the odd parts of each term:")
    print(f"odd_part({terms[0]}) * odd_part({terms[1]}) * odd_part({terms[2]}) * odd_part({terms[3]})")
    odd_parts_product = ' * '.join([str(t // (t & -t)) for t in terms]) # clever way to get odd part
    
    final_result = 1
    result_str_parts = []
    for p, exp in sorted_factors:
        final_result *= (p**exp)
        result_str_parts.append(f"{p}^{exp}")
        
    equation_lhs = " * ".join(result_str_parts)

    print(f"Highest possible order for E = {equation_lhs} = {final_result}")

solve()
<<<315>>>