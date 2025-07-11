import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of num.
    e.g., get_prime_factorization(12) -> {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp_num = num
    while d * d <= temp_num:
        while temp_num % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_num //= d
        d += 1
    if temp_num > 1:
        factors[temp_num] = factors.get(temp_num, 0) + 1
    return factors

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    - mu(n) = 1 if n is a square-free positive integer with an even number of prime factors.
    - mu(n) = -1 if n is a square-free positive integer with an odd number of prime factors.
    - mu(n) = 0 if n has a squared prime factor.
    """
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
            
    if len(factors) % 2 == 1:
        return -1
    else:
        return 1

def solve_bch_order_10():
    """
    Calculates the number of nonzero coefficients of order 10 in the BCH expansion.
    """
    n = 10
    k = 2

    # Find divisors of n
    divisors = [d for d in range(1, n + 1) if n % d == 0]
    
    total_sum = 0
    
    # Store terms for printing
    equation_parts = []
    sum_parts = []

    for d in divisors:
        mu_d = mobius(d)
        n_over_d = n // d
        term = mu_d * (k ** n_over_d)
        total_sum += term
        
        sign = "+" if term >= 0 else "-"
        # Don't show + for the first term
        if not equation_parts:
            sign = "" if term >= 0 else "-"

        equation_parts.append(f"{sign} {abs(mu_d)}*({k}^({n}//{d}))")
        sum_parts.append(f"{sign} {abs(term)}")
        
    result = total_sum // n
    
    print(f"Witt's formula for n={n} and k={k} is: N = (1/{n}) * Σ_{{d|{n}}} μ(d) * k^(n/d)")
    print("\nThe divisors of 10 are: 1, 2, 5, 10.")
    
    print("\nCalculation:")
    
    # Build the full equation string with values substituted
    equation_str = f"N = (1/{n}) * [ " + " ".join(equation_parts) + " ]"
    equation_str = equation_str.replace('[ +', '[')
    print(equation_str)

    sum_str = f"N = (1/{n}) * [ " + " ".join(sum_parts) + " ]"
    sum_str = sum_str.replace('[ +', '[')
    print(sum_str)
    
    print(f"N = (1/{n}) * [ {total_sum} ]")
    print(f"N = {result}")

    print("\nThus, the number of nonzero coefficients of order 10 is:")
    print(result)

solve_bch_order_10()
