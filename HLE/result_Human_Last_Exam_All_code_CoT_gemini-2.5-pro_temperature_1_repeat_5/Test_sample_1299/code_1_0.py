import math

def get_prime_factorization(n):
    """
    Calculates the prime factorization of a given integer n.
    Returns a dictionary where keys are prime factors and values are their exponents.
    """
    if n <= 0:
        raise ValueError("Input must be a positive integer.")
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

def print_cardinalities(l):
    """
    Calculates and prints the cardinalities of U_l and T_l based on derived formulas.
    """
    print(f"--- Calculating for l = {l} ---")
    
    # Handle the special case of l=1
    if l == 1:
        # Part A
        print("A) For l=1, there are no prime factors, so s=0.")
        print("   |U_1| = 2^s = 2^0 = 1")
        
        # Part B
        print("B) For l=1, the only triple is (m,n,lambda)=(1,1,1).")
        print("   The condition is lambda^2 < max(1^2, 2), which is 1 < 2 (True).")
        print("   |T_1| = 1")
        return

    # For l > 1
    factors = get_prime_factorization(l)
    primes = list(factors.keys())
    exponents = list(factors.values())
    s = len(primes)

    # Part A
    U_l_card = 2**s
    factor_str = ' * '.join([f'{p}^{e}' for p, e in factors.items()])
    print(f"A) Prime factorization of {l} is {factor_str}.")
    print(f"   Number of distinct prime factors, s = {s}.")
    print(f"   Formula: |U_l| = 2^s")
    print(f"   |U_{l}| = 2^{s} = {U_l_card}")

    # Part B
    prod_val = 1
    prod_str_list = []
    calc_str_terms = []

    for e in exponents:
        term = 1 + 2 * e
        prod_val *= term
        prod_str_list.append(f"(1 + 2*{e})")
        calc_str_terms.append(str(term))
    
    T_l_card = prod_val - 1
    
    print(f"B) Exponents of prime factors are e_i = {exponents}.")
    print(f"   Formula: |T_l| = (product of (1 + 2*e_i)) - 1")
    print(f"   |T_{l}| = {' * '.join(prod_str_list)} - 1")
    if s > 1:
        print(f"         = {' * '.join(calc_str_terms)} - 1")
    print(f"         = {prod_val} - 1 = {T_l_card}")

# Example with a user-provided value, e.g., l = 60
l_value = 60
print_cardinalities(l_value)

# Print the final answer in the requested format, assuming l > 1
# This part is for the final submission format, using the general symbolic expressions
s_str = "s"
e_i_str = "e_1, e_2, ..., e_s"
A_expr = f"2^s"
B_expr = f"(1+2*e_1)*(1+2*e_2)*...*(1+2*e_s) - 1"
#print("\nFinal Answer Format:")
#print(f"<<<A)[{A_expr}] B)[{B_expr}]>>>")