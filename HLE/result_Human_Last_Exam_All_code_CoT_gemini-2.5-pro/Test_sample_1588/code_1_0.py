import math

def get_prime_factorization(num):
    """
    Returns the prime factorization of a number as a dictionary.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """
    Calculates the Möbius function μ(n).
    μ(1) = 1
    μ(n) = 0 if n has a squared prime factor.
    μ(n) = (-1)^k if n is a product of k distinct primes.
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

def solve_bch_coeffs():
    """
    Calculates the number of nonzero BCH coefficients for a given order and number of generators.
    """
    n = 10  # The order
    k = 2   # The number of generators (X and Y)

    print(f"Calculating the number of nonzero BCH coefficients of order n={n} for k={k} generators.")
    print("Using Witt's formula: L_k(n) = (1/n) * sum_{d|n} (mu(d) * k^(n/d))")
    print("-" * 20)

    divisors = get_divisors(n)
    print(f"The divisors 'd' of n={n} are: {divisors}")

    sum_terms_values = []
    sum_terms_str = []
    
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term_value = mu_d * (k ** power)
        
        sum_terms_values.append(term_value)
        sum_terms_str.append(f"(mu({d}) * {k}^({n}/{d}))")

    print(f"\nL_{k}({n}) = (1/{n}) * [ {' + '.join(sum_terms_str)} ]")

    # Show the values of mu(d) and the powers
    sum_terms_calc_str = []
    for i, d in enumerate(divisors):
        mu_d = mobius(d)
        power = n // d
        sum_terms_calc_str.append(f"({mu_d} * {k}^{power})")
    
    print(f"L_{k}({n}) = (1/{n}) * [ {' + '.join(sum_terms_calc_str).replace('+ -', '- ')} ]")
    
    # Show the result of each term
    sum_terms_results_str = [str(val) for val in sum_terms_values]
    print(f"L_{k}({n}) = (1/{n}) * [ {' + '.join(sum_terms_results_str).replace('+ -', '- ')} ]")

    # Show the final sum
    total_sum = sum(sum_terms_values)
    print(f"L_{k}({n}) = (1/{n}) * [ {total_sum} ]")
    
    # Show the final result
    result = total_sum // n
    print(f"L_{k}({n}) = {result}")
    
    print("\nTherefore, the number of nonzero coefficients is " + str(result) + ".")


solve_bch_coeffs()
<<<99>>>