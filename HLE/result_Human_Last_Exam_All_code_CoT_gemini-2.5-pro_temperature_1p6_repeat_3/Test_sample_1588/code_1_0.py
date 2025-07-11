import math

def get_prime_factorization(num):
    """Finds the prime factorization of a number."""
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
    """Calculates the Mobius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def get_divisors(n):
    """Gets all positive divisors of a number."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_order_10():
    """
    Calculates the number of nonzero coefficients of order 10 in BCH expansion
    using Witt's formula.
    """
    n = 10  # The order
    k = 2   # The number of generators (X and Y)

    divisors = get_divisors(n)
    
    terms = []
    total_sum = 0

    # Build the equation string parts
    sum_str_parts = []
    calc_str_parts = []
    
    for d in divisors:
        mu_d = mobius(d)
        power_val = k**(n // d)
        term_val = mu_d * power_val
        total_sum += term_val
        
        sum_str_parts.append(f"μ({d})*({k}^({n}//{d}))")
        calc_str_parts.append(f"({mu_d})*({power_val})")

    result = total_sum // n

    # Print the steps of the calculation
    print(f"Witt's Formula for n={n}, k={k}: L_k(n) = (1/n) * Σ_{{d|n}} μ(d) * k^(n/d)")
    print(f"Divisors of {n}: {divisors}")
    print("\nCalculating the sum part of the formula:")
    equation = f"L_{k}({n}) = (1/{n}) * [ {' + '.join(sum_str_parts)} ]"
    print(equation)
    
    calculation = f"L_{k}({n}) = (1/{n}) * [ {' + '.join(calc_str_parts)} ]"
    # Replace '+ (-' with '- ' for better readability
    print(calculation.replace('+ (-', '- '))

    intermediate_sum = f"L_{k}({n}) = (1/{n}) * [ { ' - '.join(map(str, [t for t in [1024, 32, 4] if t > 0])) } + 2 ]"
    final_sum = f"L_{k}({n}) = (1/{n}) * [ {1024 - 32 - 4 + 2} ]"
    print(f"L_{k}({n}) = (1/{n}) * [ {1024} - {32} - {4} + {2} ]")
    print(final_sum)
    print(f"L_{k}({n}) = {total_sum} / {n}")
    print(f"\nThe number of nonzero coefficients is: {result}")


solve_bch_order_10()
<<<99>>>