import math

def prime_factors(n):
    """
    Returns a dictionary of prime factors of n and their multiplicities.
    e.g., prime_factors(12) returns {2: 2, 3: 1}
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
    """
    if n == 1:
        return 1
    
    factors = prime_factors(n)
    for p in factors:
        if factors[p] > 1:
            return 0  # n has a squared prime factor
            
    # n is square-free
    return (-1)**len(factors)

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_order_10():
    """
    Calculates the number of nonzero coefficients of order 10 in BCH expansion
    and prints the calculation step-by-step.
    """
    n = 10
    divisors = get_divisors(n)
    
    print(f"The number of nonzero coefficients of order {n} in the Baker-Campbell-Hausdorff expansion is given by Witt's formula:")
    print(f"L({n}) = (1/{n}) * Σ [over d|{n}] mu(d) * 2^({n}/d)\n")
    
    # Build string for the formula with symbols
    formula_symbolic_parts = []
    for d in divisors:
        formula_symbolic_parts.append(f"μ({d})*2^({n}/{d})")
    print(f"L({n}) = (1/{n}) * ( " + " + ".join(formula_symbolic_parts) + " )")
    
    # Build string for the formula with calculated values and perform calculation
    total_sum = 0
    formula_values_parts = []
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term = mu_d * (2**power)
        total_sum += term
        formula_values_parts.append(f"({mu_d})*2^{power}")
        
    print(f"L({n}) = (1/{n}) * ( " + " + ".join(formula_values_parts) + " )")

    # Build string for the detailed calculation
    calculation_parts = []
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        base_power_val = 2**power
        calculation_parts.append(f"({mu_d})*{base_power_val}")

    calc_str = " + ".join(calculation_parts).replace("+ (-", "- ")
    print(f"L({n}) = (1/{n}) * ( {calc_str} )")
    
    # Final steps
    print(f"L({n}) = (1/{n}) * ( {total_sum} )")
    result = total_sum // n
    print(f"L({n}) = {result}")
    
    print(f"\nThus, the number of nonzero coefficients of order 10 is {result}.")

if __name__ == "__main__":
    solve_bch_order_10()
<<<99>>>