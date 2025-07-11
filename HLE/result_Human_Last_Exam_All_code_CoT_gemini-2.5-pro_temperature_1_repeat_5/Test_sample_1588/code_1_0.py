import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors and their powers."""
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
    """Calculates the Möbius function μ(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_order_10():
    """
    Calculates and prints the number of nonzero coefficients of order 10
    in the Baker-Campbell-Hausdorff expansion.
    """
    n = 10  # order
    k = 2   # number of generators (X and Y)

    divisors = get_divisors(n)
    total_sum = 0
    
    final_eq_parts = []
    
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        base_val = k ** power
        term_val = mu_d * base_val
        total_sum += term_val
        final_eq_parts.append(f"({mu_d}) * {base_val}")

    result = total_sum // n
    
    final_equation = f"(1/{n}) * ({' + '.join(final_eq_parts).replace('+ (-', '- (')}) = {result}"

    print("The number of nonzero coefficients is calculated using Witt's formula for a free Lie algebra.")
    print("For order n=10 and k=2 generators, the formula is:")
    print(f"L_{n}({k}) = (1/{n}) * Σ_{{d|{n}}} μ(d) * {k}^({n}/d)")
    print(f"\nThe divisors of {n} are {divisors}.")
    print("The corresponding values for the Möbius function μ(d) are μ(1)=1, μ(2)=-1, μ(5)=-1, μ(10)=1.")
    print("\nSubstituting these values into the formula gives the final equation:")
    print(final_equation)

solve_bch_order_10()