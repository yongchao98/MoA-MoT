import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors and their exponents.
    Example: get_prime_factorization(12) returns {2: 2, 3: 1}
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

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    mu(1) = 1
    mu(n) = 0 if n has a squared prime factor.
    mu(n) = (-1)^k if n is a product of k distinct primes.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def calculate_bch_nonzero_coefficients():
    """
    Calculates the number of non-zero coefficients for a given order in the BCH expansion.
    This is determined by applying Witt's formula for a free Lie algebra.
    """
    # The order of the term in the BCH expansion
    n = 10
    # The number of generators (typically X and Y)
    k = 2

    print(f"Finding the number of non-zero coefficients of order n={n}.")
    print("This is calculated using Witt's formula for k=2 generators:")
    print(f"L_k(n) = (1/n) * Sum over d|n [ mu(d) * k^(n/d) ]\n")

    divisors = get_divisors(n)
    print(f"The divisors of n={n} are: {divisors}\n")
    
    total_sum = 0
    equation_terms = []
    
    # Calculate each term in the summation
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        k_power_term = k ** power
        term = mu_d * k_power_term
        total_sum += term
        
        # Store components for the final equation printout
        equation_terms.append(f"({mu_d} * {k}^{power})")

    result = total_sum // n

    # As requested, printing the final equation with all numbers
    print("The full calculation using these values is:")
    sum_str = " + ".join(equation_terms).replace('+ (-', '- ')
    
    # Calculate numeric values for the final display
    numeric_terms = []
    current_sum = 0
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        k_power_term = k ** power
        term = mu_d * k_power_term
        numeric_terms.append(str(term))
        current_sum += term

    sum_numeric_str = " + ".join(numeric_terms).replace('+ -', '- ')
    
    print(f"Number = (1/{n}) * [ {sum_numeric_str} ]")
    print(f"Number = (1/{n}) * [ {current_sum} ]")
    print(f"Number = {result}")
    
    print("\nTherefore, the number of nonzero coefficients of order 10 is:")
    print(result)

if __name__ == '__main__':
    calculate_bch_nonzero_coefficients()
