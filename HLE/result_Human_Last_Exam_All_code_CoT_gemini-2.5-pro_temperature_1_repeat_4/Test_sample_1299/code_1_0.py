import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors and their exponents.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
       factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def solve_cardinality(l):
    """
    Calculates the cardinality of T_l and U_l.
    """
    print(f"For l = {l}:")
    
    if l <= 0:
        print("l must be a positive integer.")
        return

    # Part A: |U_l| is stated to be equal to |T_l|
    # Part B: We calculate |T_l|
    if l == 1:
        result = 1
        print("A) |U_1| = 1")
        print("B) |T_1| = 1")
        print("Final Answer: 1")
        return

    # Case l > 1
    prime_factors = get_prime_factorization(l)
    exponents = list(prime_factors.values())
    
    # Calculate the product term
    product_term = 1
    prod_str_parts = []
    for e in exponents:
        term = 2 * e + 1
        product_term *= term
        prod_str_parts.append(f"(2*{e}+1)")
        
    result = (product_term - 1) / 2
    
    # The variable d is 0 if l is even, 1 if l is odd.
    # Our derived formula does not depend on d.
    # s is the number of distinct prime factors
    s = len(exponents)
    
    # Format the expression string
    expression_str = f"1/2 * ( (product_{{i=1}}^{s}(2*e_i+1)) - 1 )"
    
    # Format the calculation string
    calc_str = f"1/2 * (({' * '.join(prod_str_parts)}) - 1)"
    if len(prod_str_parts) == 1:
        calc_str = f"1/2 * (({prod_str_parts[0]}) - 1)"

    # Print the results as requested
    print(f"A) The cardinality |U_l| is equal to |T_l|.")
    print(f"B) The formula for |T_l| where l = p_1^e_1 * ... * p_s^e_s > 1 is:")
    print(f"   |T_l| = {expression_str}")
    print(f"   For l = {l}, the exponents are {exponents}")
    print(f"   Calculation: {calc_str} = 1/2 * ({product_term} - 1) = {int(result)}")
    print(f"Final Answer: {int(result)}")

# Example usage with a value for l
l = 12
solve_cardinality(l)

# Another example
# l = 9
# solve_cardinality(l)

# Example with l=1
# l = 1
# solve_cardinality(l)