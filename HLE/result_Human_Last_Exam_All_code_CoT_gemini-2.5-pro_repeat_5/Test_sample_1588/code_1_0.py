import math

def get_prime_factorization(num):
    """
    Returns the prime factorization of a number as a dictionary.
    e.g., get_prime_factorization(12) -> {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    - mu(1) = 1
    - mu(n) = 0 if n has a squared prime factor.
    - mu(n) = (-1)^k if n is a product of k distinct primes.
    """
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
            
    return (-1)**len(factors)

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_coeffs(n, k):
    """
    Calculates the number of nonzero coefficients for BCH expansion of order n
    with k generators using Witt's formula and prints the steps.
    """
    divisors = get_divisors(n)
    
    sum_terms = []
    sum_values = []
    
    total_sum_val = 0
    
    for d in divisors:
        mu_d = mobius(d)
        power_val = n // d
        term_val = mu_d * (k**power_val)
        
        # Build string for each term in the sum
        term_str = f"({mu_d} * {k}^{power_val})"
        sum_terms.append(term_str)
        
        # Keep track of calculated values
        sum_values.append(str(term_val))
        total_sum_val += term_val

    # Print the formula with symbols
    formula_str = f"L({n}, {k}) = (1/{n}) * ( " + " + ".join(sum_terms) + " )"
    print(formula_str)

    # Print the formula with intermediate calculated values
    values_str = f"L({n}, {k}) = (1/{n}) * ( " + " + ".join(sum_values).replace('+ -', '- ') + " )"
    print(values_str)

    # Print the formula with the sum evaluated
    sum_result_str = f"L({n}, {k}) = (1/{n}) * {total_sum_val}"
    print(sum_result_str)
    
    # Calculate and print the final result
    final_result = total_sum_val // n
    final_result_str = f"L({n}, {k}) = {final_result}"
    print(final_result_str)
    
    print(f"\nThe number of nonzero coefficients of order {n} is: {final_result}")
    
    return final_result

# Parameters for the problem
order_n = 10
generators_k = 2

# Solve the problem
result = solve_bch_coeffs(order_n, generators_k)