import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors of num and their powers."""
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

def moebius(m):
    """Calculates the Moebius function mu(m)."""
    if m == 1:
        return 1
    factors = get_prime_factorization(m)
    for p in factors:
        if factors[p] > 1:
            return 0  # mu(m) = 0 if m has a squared prime factor
    return (-1) ** len(factors)

def get_divisors(m):
    """Returns a sorted list of all divisors of m."""
    divs = set()
    for i in range(1, int(math.sqrt(m)) + 1):
        if m % i == 0:
            divs.add(i)
            divs.add(m // i)
    return sorted(list(divs))

def solve_bch_count():
    """
    Calculates the number of nonzero BCH coefficients for a given order
    using Witt's formula.
    """
    n = 10
    k = 2 # Number of generators (X, Y)

    print(f"The number of nonzero coefficients of order {n} in the BCH expansion is given by Witt's formula:")
    print(f"N(n) = (1/n) * sum_{{d|n}} mu(n/d) * k^d")
    print(f"Here, n = {n} and k = {k}.")
    print("-" * 20)

    divs = get_divisors(n)
    print(f"The divisors 'd' of {n} are: {divs}")

    total_sum = 0
    sum_expr_parts = []
    
    # Building the expression string for the numerator sum
    numerator_expr_str = "Sum = " + " + ".join([f"mu({n}/{d}) * {k}^{d}" for d in divs])
    print(f"\nCalculating the sum part of the formula:")
    print(numerator_expr_str)

    # Calculate values and build the detailed expression
    mu_values = [moebius(n // d) for d in divs]
    k_pow_d_values = [k**d for d in divs]
    term_values = [mu * kd for mu, kd in zip(mu_values, k_pow_d_values)]

    values_str = "    = " + " + ".join([f"({mu}) * {kd}" for mu, kd in zip(mu_values, k_pow_d_values)])
    print(values_str)
    
    calculation_str = "    = " + " + ".join([str(val) for val in term_values]).replace('+ -', '- ')
    total_sum = sum(term_values)
    print(f"{calculation_str}")
    print(f"    = {total_sum}")

    # Final calculation
    result = total_sum // n
    print("\nFinally, divide the sum by n:")
    print(f"Result = (1/{n}) * {total_sum}")
    print(f"Result = {result}")

    print("\nThus, the number of nonzero coefficients is:")
    print(result)

solve_bch_count()
<<<99>>>