import math

def get_prime_factorization(num):
    """Returns the prime factorization of a number as a dictionary."""
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
    Calculates the Möbius function mu(n).
    - mu(1) = 1
    - mu(n) = 0 if n has a squared prime factor.
    - mu(n) = (-1)^k if n is a product of k distinct primes.
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

def calculate_bch_coefficients_count(n, k):
    """
    Calculates the number of nonzero BCH coefficients of order n for k generators
    using Witt's formula.
    """
    divisors = get_divisors(n)
    
    numerator_sum = 0
    numerator_terms = []
    
    print(f"Calculating the number of nonzero coefficients for order n={n} and k={k} generators.")
    print("The formula is: (1/n) * Sum_{d|n} [mu(d) * k^(n/d)]")
    print(f"\nThe divisors of {n} are: {divisors}")
    print("\nCalculating each term in the sum:")

    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term = mu_d * (k**power)
        
        numerator_sum += term
        numerator_terms.append(term)
        
        print(f"  d={d}: mu({d}) * {k}^({n}/{d}) = {mu_d} * {k}^{power} = {term}")

    result = numerator_sum // n
    
    # Construct the final equation string
    sum_str = " + ".join(str(t) for t in numerator_terms).replace("+ -", "- ")
    
    print("\nFinal Equation:")
    print(f"({sum_str}) / {n} = {numerator_sum} / {n} = {result}")

    print(f"\nThe number of nonzero coefficients of order {n} is {result}.")


# Parameters for the problem
order_n = 10
generators_k = 2

calculate_bch_coefficients_count(order_n, generators_k)