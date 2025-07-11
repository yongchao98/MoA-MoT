import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors and their exponents.
    e.g., get_prime_factorization(12) -> {2: 2, 3: 1}
    """
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

def calculate_cardinality_T_ell(l):
    """
    Calculates the cardinality of T_l based on the derived formula.
    |T_l| = 1 if l=1
    |T_l| = (product of (1+2e_i)) - 1 if l > 1
    where l = p1^e1 * p2^e2 * ...
    """
    if not isinstance(l, int) or l <= 0:
        raise ValueError("l must be a positive integer.")
    
    if l == 1:
        return 1, "1", "1"
    
    prime_factors = get_prime_factorization(l)
    exponents = list(prime_factors.values())
    
    product_term = 1
    prod_str = ""
    for e in exponents:
        term = 1 + 2 * e
        product_term *= term
        prod_str += f"(1+2*{e})"
    
    result = product_term - 1
    
    # Building the formula string
    # Remove single outermost parens if they exist, not needed for single factor
    if len(exponents) == 1:
        prod_str = prod_str.strip('()')
        
    formula_str = f"{prod_str} - 1"
    
    # Calculate intermediate steps for display
    terms_calc = [1 + 2 * e for e in exponents]
    calc_str_parts = []
    if len(terms_calc) > 1:
        for t in terms_calc:
             calc_str_parts.append(str(t))
        calc_str = f"({'*'.join(calc_str_parts)}) - 1 = {product_term} - 1"
    else:
        calc_str = f"{terms_calc[0]} - 1"
    
    final_calc_str = f"{formula_str} = {calc_str} = {result}" if len(terms_calc)>0 else "1"


    return result, formula_str, final_calc_str

def main():
    """
    Main function to solve the task.
    You can change the value of l to test with other numbers.
    """
    # Example value for l
    l = 12 

    # Prime factorization of l=12 is 2^2 * 3^1
    # The exponents are e1=2, e2=1.
    
    cardinality, formula, final_calc_str = calculate_cardinality_T_ell(l)

    print(f"For l = {l}:")
    print(f"The cardinality is |T_{l}| = {cardinality}.")
    print(f"The calculation from the formula is: {final_calc_str}")

if __name__ == "__main__":
    main()
