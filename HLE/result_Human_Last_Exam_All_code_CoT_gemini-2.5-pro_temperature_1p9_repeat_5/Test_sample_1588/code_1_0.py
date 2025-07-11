import math

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    mu(1) = 1
    mu(n) = 0 if n has a squared prime factor.
    mu(n) = (-1)^k if n is a product of k distinct primes.
    """
    if n == 1:
        return 1
    
    prime_factors = 0
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            prime_factors += 1
            temp_n //= d
            if temp_n % d == 0:  # Check for squared factor
                return 0
        d += 1
    if temp_n > 1: # Remaining factor is prime
        prime_factors += 1
        
    return -1 if prime_factors % 2 != 0 else 1

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def calculate_bch_coeffs():
    """
    Calculates the number of nonzero BCH coefficients for n=10, k=2
    using Witt's formula and prints the calculation steps.
    """
    n = 10  # Order of the expansion
    k = 2   # Number of generators

    divisors = get_divisors(n)
    
    total_sum = 0
    
    # Store parts of the equation for printing
    formula_terms = []
    evaluated_terms = []
    numeric_values = []

    for d in divisors:
        mu_d = mobius(d)
        exponent = n // d
        term_value = mu_d * (k ** exponent)
        total_sum += term_value

        # Build strings for each step of the equation
        mu_str = f"({mu_d})" if mu_d < 0 else str(mu_d)
        formula_terms.append(f"{mu_str} * {k}^{exponent}")
        evaluated_terms.append(f"{mu_str} * {k**exponent}")
        numeric_values.append(term_value)
    
    final_result = total_sum // n

    # Assemble the calculation string (e.g., "1024 - 32 - 4 + 2")
    sum_str = ""
    for val in numeric_values:
        if val >= 0 and sum_str:
            sum_str += f" + {val}"
        else:
            sum_str += f" {val}"
    sum_str = sum_str.strip()

    # Print the step-by-step calculation
    print(f"The number of nonzero coefficients of order {n}, L_{n}({k}), is calculated using Witt's formula:")
    print(f"L_{n}({k}) = (1/{n}) * ( " + " + ".join(formula_terms) + " )")
    print(f"L_{n}({k}) = (1/{n}) * ( " + " + ".join(evaluated_terms) + " )")
    print(f"L_{n}({k}) = (1/{n}) * ( {sum_str} )")
    print(f"L_{n}({k}) = (1/{n}) * ( {total_sum} )")
    print(f"L_{n}({k}) = {final_result}")

calculate_bch_coeffs()