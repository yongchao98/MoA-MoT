import math

def get_prime_factorization_exponents(n):
    """
    Calculates the exponents of the prime factorization of n.
    For n = p1^e1 * p2^e2 * ... * ps^es, returns [e1, e2, ..., es].
    For n = 1, returns an empty list.
    """
    if n == 1:
        return []
    
    exponents = []
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            count = 0
            while temp_n % d == 0:
                count += 1
                temp_n //= d
            exponents.append(count)
        d += 1
    if temp_n > 1:
        exponents.append(1)
    return exponents

def calculate_T_ell(l):
    """
    Calculates the cardinality of T_l based on the derived piecewise formula.
    """
    if not isinstance(l, int) or l <= 0:
        raise ValueError("l must be a positive integer.")

    if l == 1:
        return 1
    else:
        # For l > 1, the formula is (product of (2*e_i + 1)) - 1
        exponents = get_prime_factorization_exponents(l)
        
        # This is tau(l^2)
        prod_val = 1
        for e in exponents:
            prod_val *= (2 * e + 1)
            
        return prod_val - 1

def main():
    # Example for l = p1^e1 * ...
    # This example demonstrates the calculation for l=6.
    # l = 6 = 2^1 * 3^1 => e1=1, e2=1
    l = 6
    print(f"For l = {l}:")
    # According to the formula for l > 1:
    # It is ( (2*e1 + 1) * (2*e2 + 1) * ... ) - 1
    # For l=6, e1=1, e2=1. The result is (2*1+1)*(2*1+1) - 1 = 3*3 - 1 = 8.
    
    # Let's verify with the code.
    result = calculate_T_ell(l)
    
    # We construct the equation string for clarity
    exponents = get_prime_factorization_exponents(l)
    tau_l_squared_terms = [f"(2*{e}+1)" for e in exponents]
    tau_l_squared_expr = "*".join(tau_l_squared_terms)
    tau_l_squared_val = 1
    for e in exponents:
        tau_l_squared_val *= (2 * e + 1)
    
    if l == 1:
        print(f"|T_1| = 1")
    else:
        print(f"|T_{l}| = ({tau_l_squared_expr}) - 1 = {tau_l_squared_val} - 1 = {result}")

    # Another example for l=1
    l = 1
    result_1 = calculate_T_ell(l)
    print(f"\nFor l = {l}:")
    print(f"|T_1| = {result_1}")
    
    # Another example for l=9
    l = 9 # 3^2, so e1=2
    result_9 = calculate_T_ell(l)
    print(f"\nFor l = {l}:")
    exponents_9 = get_prime_factorization_exponents(l)
    tau_l_squared_terms_9 = [f"(2*{e}+1)" for e in exponents_9]
    tau_l_squared_expr_9 = "*".join(tau_l_squared_terms_9)
    tau_l_squared_val_9 = 1
    for e in exponents_9:
        tau_l_squared_val_9 *= (2 * e + 1)
    print(f"|T_{l}| = ({tau_l_squared_expr_9}) - 1 = {tau_l_squared_val_9} - 1 = {result_9}")

if __name__ == "__main__":
    main()
