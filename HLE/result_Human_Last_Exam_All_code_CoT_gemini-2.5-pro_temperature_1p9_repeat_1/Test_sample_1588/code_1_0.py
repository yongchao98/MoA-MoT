import math

def get_prime_factorization(num):
    """
    Returns a list of prime factors of a given number.
    Example: get_prime_factorization(12) -> [2, 2, 3]
    """
    factors = []
    d = 2
    temp_num = num
    while d * d <= temp_num:
        while temp_num % d == 0:
            factors.append(d)
            temp_num //= d
        d += 1
    if temp_num > 1:
        factors.append(temp_num)
    return factors

def mobius(n):
    """
    Calculates the Möbius function mu(n).
    - mu(n) = 1 if n is a square-free positive integer with an even number of prime factors.
    - mu(n) = -1 if n is a square-free positive integer with an odd number of prime factors.
    - mu(n) = 0 if n has a squared prime factor.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    if n == 1:
        return 1
    
    prime_factors = get_prime_factorization(n)
    
    # Check for squared prime factors
    if len(prime_factors) != len(set(prime_factors)):
        return 0
    else:
        # It's square-free
        return (-1)**len(prime_factors)

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of a number n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def calculate_bch_nonzero_coeffs(n, k):
    """
    Calculates the number of nonzero coefficients for a given order n
    in the BCH expansion with k generators using Witt's formula and
    prints the calculation steps.
    """
    divisors = get_divisors(n)
    total_sum = 0
    
    print("The number of nonzero coefficients is calculated using Witt's formula:")
    print(f"  Number = (1/n) * Σ_{{d|n}} μ(n/d) * k^d\n")
    print(f"For this problem, n = {n} (order 10) and k = 2 (generators X and Y).")
    
    divisors = get_divisors(n)
    print(f"\nThe divisors 'd' of n={n} are: {divisors}")
    
    sum_terms_values = []
    
    print("\nCalculating each term 'μ(n/d) * k^d' in the sum:")
    for d in divisors:
        n_over_d = n // d
        mu_val = mobius(n_over_d)
        term_val = mu_val * (k**d)
        total_sum += term_val
        sum_terms_values.append(term_val)
        
        print(f"  For d={d:2d}: μ(n/d) = μ({n}/{d}) = μ({n_over_d}) = {mu_val:2d}. The term is {mu_val:2d} * {k}^{d} = {term_val}")
    
    # Construct the final equation string for printing
    equation_sum_str = ' + '.join(map(str, sum_terms_values)).replace('+ -', '- ')
    
    final_result = total_sum // n
    
    print("\nFinal calculation:")
    print(f"Number = (1/{n}) * ({equation_sum_str})")
    print(f"Number = (1/{n}) * ({total_sum})")
    print(f"Number = {final_result}")

if __name__ == "__main__":
    # Define the order and number of generators
    order = 10
    generators = 2
    
    calculate_bch_nonzero_coeffs(order, generators)