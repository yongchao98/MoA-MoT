import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary of {prime: exponent}.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def solve_task():
    """
    Calculates the number of power subgroups in the generalized quaternion group of size 128.
    """
    # The exponent of the generalized quaternion group Q_128 is 64.
    # The number of power subgroups is the number of divisors of the exponent.
    exponent = 64

    print(f"The exponent of the generalized quaternion group of size 128 is {exponent}.")
    print(f"The number of power subgroups is equal to the number of divisors of {exponent}.")
    
    # Find the prime factorization of the exponent
    factors = get_prime_factorization(exponent)
    
    # Build the explanation strings
    factor_str_parts = []
    for prime, power in factors.items():
        factor_str_parts.append(f"{prime}^{power}")
    factor_str = " * ".join(factor_str_parts)
    
    print(f"\nTo find the number of divisors, we first get the prime factorization of {exponent}:")
    print(f"{exponent} = {factor_str}")

    # Use the formula for the number of divisors: (e1+1)*(e2+1)*...
    num_divisors = 1
    formula_parts = []
    result_parts = []
    for prime, power in factors.items():
        num_divisors *= (power + 1)
        formula_parts.append(f"({power} + 1)")
        result_parts.append(str(power + 1))

    formula_str = " * ".join(formula_parts)
    result_str = " * ".join(result_parts)

    print("\nThe number of divisors is calculated from the exponents in the factorization:")
    if len(result_parts) > 1:
        print(f"Number of divisors = {formula_str} = {result_str} = {num_divisors}")
    else:
        print(f"Number of divisors = {formula_str} = {num_divisors}")

solve_task()