import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of num.
    e.g., get_prime_factorization(12) -> {2: 2, 3: 1}
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
    Computes the Mobius function μ(n).
    μ(n) = 1 if n is a square-free positive integer with an even number of prime factors.
    μ(n) = -1 if n is a square-free positive integer with an odd number of prime factors.
    μ(n) = 0 if n has a squared prime factor.
    """
    if n == 1:
        return 1
    
    prime_factors = get_prime_factorization(n)
    # If n has a squared prime factor, μ(n) = 0
    for p in prime_factors:
        if prime_factors[p] > 1:
            return 0
    # Otherwise, (-1)^(number of distinct prime factors)
    return (-1) ** len(prime_factors)

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_coefficients():
    """
    Calculates the number of nonzero BCH coefficients for a given order and number of generators.
    """
    n = 10  # The order of the expansion
    k = 2   # The number of generators (e.g., X and Y)

    print(f"Calculating the number of nonzero coefficients of order n={n} in the BCH expansion for k={k} generators.")
    print("Using Witt's formula: L_n(k) = (1/n) * Σ_{d|n} μ(d) * k^(n/d)\n")
    
    divisors_of_n = get_divisors(n)
    print(f"The divisors of n={n} are: {divisors_of_n}")
    
    total_sum = 0
    sum_parts = []
    
    for d in divisors_of_n:
        mu_d = mobius(d)
        exponent = n // d
        term_value = mu_d * (k ** exponent)
        total_sum += term_value
        sum_parts.append(f"(μ({d})) * {k}^({n}/{d})")

    # Displaying the symbolic formula
    full_formula = f"L_{n}({k}) = (1/{n}) * [ {' + '.join(sum_parts)} ]"
    print(f"\nFormula breakdown:\n{full_formula}")

    # Displaying the intermediate values
    value_parts = []
    for d in divisors_of_n:
        mu_d = mobius(d)
        exponent = n // d
        base_val = k ** exponent
        value_parts.append(f"({mu_d})*({base_val})")

    full_values = f"L_{n}({k}) = (1/{n}) * [ {' + '.join(value_parts)} ]"
    print(f"\nPlugging in the numbers:\n{full_values}")

    # Displaying the final calculation
    final_calc_str = f"L_{n}({k}) = (1/{n}) * [ {total_sum} ]"
    print(f"\nCalculating the sum:\n{final_calc_str}")

    result = total_sum // n
    print(f"\nFinal result: L_{n}({k}) = {result}")

solve_bch_coefficients()
<<<99>>>