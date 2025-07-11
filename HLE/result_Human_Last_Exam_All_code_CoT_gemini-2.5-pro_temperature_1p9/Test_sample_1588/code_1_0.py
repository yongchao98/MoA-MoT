import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors of n.
    Example: get_prime_factorization(12) returns {2: 2, 3: 1}
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

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    mu(n) = 1 if n is a square-free positive integer with an even number of prime factors.
    mu(n) = -1 if n is a square-free positive integer with an odd number of prime factors.
    mu(n) = 0 if n has a squared prime factor.
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
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def main():
    """
    Calculates and prints the number of nonzero BCH coefficients of a given order.
    """
    n = 10  # Order of the expansion
    k = 2   # Number of generators (X and Y)

    print(f"Calculating the number of nonzero BCH coefficients of order {n}.")
    print("This is given by Witt's formula: l(n) = (1/n) * sum_{d|n}(mu(n/d) * k^d)")
    
    divisors = get_divisors(n)
    print(f"\nThe divisors 'd' of n={n} are: {divisors}")

    total_sum = 0
    sum_terms = []
    
    print("\nCalculating the terms mu(n/d) * k^d for each divisor d:")
    for d in divisors:
        n_over_d = n // d
        mu_val = mobius(n_over_d)
        term_val = mu_val * (k**d)
        total_sum += term_val
        sum_terms.append(str(term_val))
        print(f"  d={d:2}: mu({n}/{d})*k^d = mu({n_over_d})*{k}^{d} = {mu_val:2} * {k**d:<4} = {term_val}")
    
    # Ensure integer division
    result = total_sum // n

    sum_str = " + ".join(sum_terms).replace("+ -", "- ")

    print("\nPutting the terms together into the final equation:")
    print(f"l({n}) = (1/{n}) * ({sum_str})")
    print(f"l({n}) = (1/{n}) * ({total_sum})")
    print(f"l({n}) = {result}")

if __name__ == "__main__":
    main()