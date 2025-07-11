import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors of num."""
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
    """Calculates the Mobius function mu(n)."""
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

def solve_bch_coeffs():
    """
    Calculates the number of nonzero coefficients of order n in the BCH expansion
    for k generators using Witt's formula.
    """
    n = 10
    k = 2

    print(f"Calculating the number of nonzero BCH coefficients of order n={n} for k={k} generators.")
    print("Using Witt's formula: L_k(n) = (1/n) * sum_{d|n} (mu(d) * k^(n/d))")
    
    divisors = get_divisors(n)
    print(f"The divisors of {n} are: {divisors}")

    total_sum = 0
    equation_parts = []

    for d in divisors:
        mu_d = mobius(d)
        power_val = k**(n // d)
        term = mu_d * power_val
        total_sum += term
        
        # Build the string for each term in the sum
        part_str = f"({mu_d} * {k}^{n//d})"
        equation_parts.append(part_str)

    result = total_sum // n

    # Build the full equation string
    sum_str = " + ".join(equation_parts)
    print("\nStep-by-step calculation:")
    print(f"L_{k}({n}) = (1/{n}) * ( {sum_str} )")
    
    # Show the evaluated terms
    evaluated_parts = []
    for d in divisors:
        mu_d = mobius(d)
        power_val = k**(n // d)
        term = mu_d * power_val
        evaluated_parts.append(str(term))
    
    evaluated_sum_str = " + ".join(evaluated_parts).replace("+ -", "- ")
    print(f"L_{k}({n}) = (1/{n}) * ( {evaluated_sum_str} )")
    print(f"L_{k}({n}) = (1/{n}) * ( {total_sum} )")
    print(f"L_{k}({n}) = {result}")
    
    print("\nFinal Answer:")
    print(result)


solve_bch_coeffs()
<<<99>>>