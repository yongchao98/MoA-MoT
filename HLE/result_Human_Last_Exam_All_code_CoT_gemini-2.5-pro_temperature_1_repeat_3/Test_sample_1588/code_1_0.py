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

def moebius(n):
    """Calculates the Moebius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a sorted list of positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def calculate_bch_term_count(n, k):
    """
    Calculates the number of non-zero coefficients of order n for k generators
    using Witt's formula.
    """
    divisors = get_divisors(n)
    
    total_sum_val = 0
    
    # For printing the formula
    formula_terms = []
    value_terms = []
    
    for d in divisors:
        mu_d = moebius(d)
        power_val = k**(n//d)
        term_val = mu_d * power_val
        total_sum_val += term_val
        
        formula_terms.append(f"({mu_d}) * {k}**({n//d})")
        value_terms.append(f"{term_val}")

    result = total_sum_val // n
    
    print("The number of non-zero coefficients is given by Witt's formula:")
    print(f"L_{n}({k}) = (1/{n}) * Σ_{{d|{n}}} μ(d) * {k}^({n}/d)\n")
    
    print("Calculation steps:")
    # Print the expanded formula with all numbers
    formula_str = " + ".join(formula_terms).replace('+ -', '- ')
    print(f"L_{n}({k}) = (1/{n}) * ( {formula_str} )")
    
    # Print the formula with evaluated powers
    value_str = " + ".join(value_terms).replace('+ -', '- ')
    print(f"L_{n}({k}) = (1/{n}) * ( {value_str} )")

    # Print the sum
    print(f"L_{n}({k}) = (1/{n}) * {total_sum_val}")
    
    # Print the final result
    print(f"L_{n}({k}) = {total_sum_val // n}")
    print("\n----------------------------------------------------")
    print(f"The number of nonzero coefficients of order {n} is: {result}")
    
# For the specific problem: n=10, k=2
n = 10
k = 2
calculate_bch_term_count(n, k)