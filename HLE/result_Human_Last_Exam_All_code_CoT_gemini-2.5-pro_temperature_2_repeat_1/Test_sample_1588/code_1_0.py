import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors and their powers."""
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
            divs.add(n//i)
    return sorted(list(divs))

def calculate_bch_term_count(n, k):
    """
    Calculates the number of nonzero BCH coefficients of order n
    using Witt's formula.
    """
    divisors = get_divisors(n)
    
    total_sum = 0
    
    print(f"The number of nonzero coefficients of order {n} is given by Witt's formula:")
    print(f"L_{k}({n}) = (1/{n}) * Σ_{{d|{n}}} μ(d) * {k}^({n}/d)\n")
    
    equation_parts = []
    calculation_parts = []
    
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term = mu_d * (k**power)
        total_sum += term
        
        # Build strings for pretty printing the equation
        mu_str = f"μ({d})"
        power_str = f"{k}^({n}/{d})"
        
        mu_val_str = f"({mu_d})" if mu_d < 0 else str(mu_d)
        val_str = f"{k}^{power}"
        
        equation_parts.append(f"{mu_str}*{power_str}")
        calculation_parts.append(f"{mu_val_str}*{val_str}")
        
    print("Substituting the divisors of 10 (1, 2, 5, 10):")
    print(f"L_2(10) = (1/10) * [ {' + '.join(equation_parts)} ]")
    print(f"L_2(10) = (1/10) * [ {' + '.join(calculation_parts)} ]")

    calc_terms = []
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        k_pow_val = k ** power
        term = mu_d * k_pow_val
        mu_val_str = f"({mu_d})" if mu_d < 0 else str(mu_d)
        calc_terms.append(f"{mu_val_str}*{k_pow_val}")

    print(f"L_2(10) = (1/10) * [ {' + '.join(calc_terms).replace('+ -', '- ')} ]")

    sum_terms = []
    for d in divisors:
        term = mobius(d) * (k**(n//d))
        sum_terms.append(str(term))

    print(f"L_2(10) = (1/10) * [ {' + '.join(sum_terms).replace('+ -', '- ')} ]")

    print(f"L_2(10) = (1/10) * ({total_sum})")
    
    result = total_sum // n
    print(f"\nFinal Answer: {result}")
    
    return result

# Parameters for the problem
order_n = 10
generators_k = 2

# Execute the calculation
final_answer = calculate_bch_term_count(order_n, generators_k)
print(f"<<<{final_answer}>>>")