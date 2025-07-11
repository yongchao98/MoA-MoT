import math

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def get_prime_factorization(n):
    """
    Returns a dictionary of the prime factors of n.
    The keys are the prime factors and the values are their powers.
    """
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

def solve_bch_coefficient_count():
    """
    Calculates and prints the number of nonzero BCH coefficients
    for a given order n and number of generators k.
    """
    n = 10  # Order of the expansion
    k = 2   # Number of generators (e.g., X and Y)

    print("The number of nonzero coefficients is given by Witt's formula:")
    print(f"L_k(n) = (1/n) * sum over d|n ( mu(d) * k^(n/d) )")
    print("\nFor order n=10 and k=2 generators:")
    
    divisors = get_divisors(n)
    total_sum = 0
    sum_terms = []
    value_terms = []

    for d in divisors:
        mu_d = mobius(d)
        power_val = n // d
        term_val = mu_d * (k**power_val)
        total_sum += term_val

        # Formatting for display
        mu_str = f"({mu_d})" if mu_d < 0 else str(mu_d)
        sum_terms.append(f"{mu_str} * {k}^{power_val}")
        value_terms.append(f"{mu_str} * {k**power_val}")

    result = total_sum // n

    # Print the step-by-step calculation
    print(f"L_2(10) = (1/{n}) * [ {' + '.join(sum_terms).replace('+ -', '- ')} ]")
    print(f"L_2(10) = (1/{n}) * [ {' + '.join(value_terms).replace('+ (', '- ').replace(')', '')} ]")
    print(f"L_2(10) = (1/{n}) * [ {1024} - {32} - {4} + {2} ]")
    print(f"L_2(10) = (1/{n}) * [ {total_sum} ]")
    print(f"L_2(10) = {result}")
    
    print(f"\nTherefore, the number of nonzero coefficients of order 10 is {result}.")

if __name__ == '__main__':
    solve_bch_coefficient_count()