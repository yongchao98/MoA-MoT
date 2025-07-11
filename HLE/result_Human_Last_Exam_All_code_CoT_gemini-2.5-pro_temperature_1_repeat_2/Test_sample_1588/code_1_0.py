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

def mobius_mu(n):
    """Calculates the Mobius function μ(n)."""
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
    Calculates the number of nonzero coefficients of order n in the BCH expansion.
    """
    n = 10
    k = 2

    divisors = get_divisors(n)
    total_sum = 0
    
    # For pretty printing the formula
    formula_terms = []
    value_terms = []

    for d in divisors:
        mu_d = mobius_mu(d)
        power = n // d
        term = mu_d * (k**power)
        total_sum += term
        
        formula_terms.append(f"μ({d})*({k}^({n}//{d}))")
        value_terms.append(f"({mu_d})*({k**power})")

    result = total_sum // n

    # Print the detailed calculation
    print(f"The number of nonzero coefficients is given by Witt's formula L_k(n) for k=2 and n=10.")
    print(f"L_2(10) = (1/{n}) * Σ_{{d|{n}}} μ(d) * k^(n/d)")
    print("\nCalculating each term for divisors of 10: {", ", ".join(map(str, divisors)), "}")
    
    formula_str = f"L_2(10) = (1/{n}) * [ " + " + ".join(formula_terms) + " ]"
    print(formula_str)

    value_str = f"L_2(10) = (1/{n}) * [ " + " + ".join(value_terms) + " ]"
    print(value_str.replace("+ -", "- "))

    calc_str = f"L_2(10) = (1/{n}) * [ {total_sum} ]"
    print(calc_str)
    
    final_str = f"L_2(10) = {result}"
    print(final_str)
    
    print(f"\nThe number of nonzero coefficients of order 10 is {result}.")

solve_bch_coeffs()
<<<99>>>