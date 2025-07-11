import math

def get_prime_factorization(n):
    """
    Returns a dictionary of the prime factorization of n.
    The keys are the prime factors and the values are their exponents.
    """
    factors = {}
    # Check for factor 2
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n = n // 2
    # Check for odd factors
    d = 3
    while d * d <= n:
        while (n % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            n = n // d
        d = d + 2
    if n > 1:
       factors[n] = factors.get(n, 0) + 1
    return factors

def calculate_cardinality_T_l(l):
    """
    Calculates the cardinality of the set T_l and prints the steps.
    """
    print(f"Calculating the cardinality for l = {l}")
    
    if not isinstance(l, int) or l <= 0:
        print("Error: l must be a positive integer.")
        return

    # Handle the special case for l = 1
    if l == 1:
        # For l=1, the condition is λ^2 < max(1^2, 2) => 1 < 2, which is true.
        # The only triple is (1, 1, 1). So |T_1| = 1.
        print("For l = 1, the cardinality |T_1| is 1.")
        return

    # For l > 1, the formula is |T_l| = (Π(2*e_i + 1)) - 1
    # where e_i are the exponents in the prime factorization of l.
    
    factors = get_prime_factorization(l)
    exponents = list(factors.values())
    
    print(f"The prime factorization of {l} is: {' * '.join([f'{p}^{e}' for p, e in factors.items()])}")
    print(f"The exponents e_i are: {', '.join(map(str, exponents))}")
    
    # Build the equation string step-by-step
    calc_strs = []
    terms = []
    for e in exponents:
        calc_strs.append(f"(2*{e}+1)")
        terms.append(2*e + 1)
    
    product = 1
    for term in terms:
        product *= term
        
    result = product - 1
    
    # Fulfilling the requirement to "output each number in the final equation"
    print("\nFinal Equation:")
    term_values_str = " * ".join(map(str, terms))
    calc_str = " * ".join(calc_strs)
    
    final_equation = f"|T_{l}| = {calc_str} - 1 = {term_values_str} - 1 = {product} - 1 = {result}"
    print(final_equation)
    
    print(f"\nThe cardinality |T_{l}| is: {result}")

# Example demonstrating the calculation for a given l. Let's use l = 60.
# l = 60 = 2^2 * 3^1 * 5^1. The exponents are e1=2, e2=1, e3=1.
# |T_60| = (2*2+1)*(2*1+1)*(2*1+1) - 1 = 5 * 3 * 3 - 1 = 45 - 1 = 44.
calculate_cardinality_T_l(60)