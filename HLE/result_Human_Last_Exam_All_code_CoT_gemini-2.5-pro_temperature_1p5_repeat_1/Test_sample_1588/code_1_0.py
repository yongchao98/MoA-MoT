import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors of num and their powers."""
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
    """
    Calculates the Mobius function mu(n).
    - mu(1) = 1
    - mu(n) = 0 if n has a squared prime factor.
    - mu(n) = (-1)^k if n is a product of k distinct primes.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    if n == 1:
        return 1
    
    prime_factors = get_prime_factorization(n)
    
    # Check for squared prime factors
    for p in prime_factors:
        if prime_factors[p] > 1:
            return 0
            
    # If no squared prime factors, mu(n) = (-1)^k where k is the number of distinct prime factors
    return (-1) ** len(prime_factors)

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def calculate_bch_nonzero_coeffs(n, k):
    """
    Calculates the number of nonzero coefficients of order n in the BCH expansion
    for k generators using Witt's formula, and prints the detailed calculation.
    """
    divisors = get_divisors(n)
    
    calculated_terms = []
    total_sum = 0
    
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        k_power_val = k ** power
        term_val = mu_d * k_power_val
        calculated_terms.append({
            'd': d, 
            'mu_d': mu_d, 
            'power': power, 
            'k_power_val': k_power_val, 
            'term_val': term_val
        })
        total_sum += term_val
        
    result = total_sum // n
    
    # --- Printing the explanation and calculation ---
    print(f"To find the number of nonzero coefficients of order {n} in the Baker-Campbell-Hausdorff expansion,")
    print("we use Witt's formula for the dimension of the n-th graded part of a free Lie algebra on k generators.")
    print("\nWitt's formula: L_k(n) = (1/n) * sum over d|n of (mu(d) * k^(n/d))")
    print(f"\nFor this problem, the order is n = {n}, and there are k = 2 generators (X and Y).")
    
    divisors_str = ", ".join(map(str, divisors))
    print(f"The positive divisors 'd' of n={n} are: {divisors_str}.")

    # Build and print the formula expression with symbolic values
    expr_parts = [f"mu({term['d']}) * {k}^({n}/{term['d']})" for term in calculated_terms]
    print(f"\nL_{k}({n}) = (1/{n}) * ( " + " + ".join(expr_parts) + " )")

    # Build and print the expression with numeric values substituted
    calc_parts_str = []
    for term in calculated_terms:
        # Show mu(d) * k^power
        calc_parts_str.append(f"({term['mu_d']}) * {k}^{term['power']}")
    print(f"L_{k}({n}) = (1/{n}) * ( " + " + ".join(calc_parts_str) + " )")
    
    calc_parts_2_str = []
    for term in calculated_terms:
        # Show mu(d) * result_of_power
        calc_parts_2_str.append(f"({term['mu_d']}) * {term['k_power_val']}")
    print(f"L_{k}({n}) = (1/{n}) * ( " + " + ".join(calc_parts_2_str) + " )")

    # Build and print the final sum expression
    final_sum_parts = []
    for i, term in enumerate(calculated_terms):
        val = term['term_val']
        if i > 0 and val >= 0:
            final_sum_parts.append(f"+ {val}")
        elif val < 0:
            final_sum_parts.append(f"- {abs(val)}")
        else:
            final_sum_parts.append(str(val))
    print(f"L_{k}({n}) = (1/{n}) * ( {' '.join(final_sum_parts)} )")
    
    # Print the sum and final result
    print(f"L_{k}({n}) = (1/{n}) * ( {total_sum} )")
    print(f"L_{k}({n}) = {total_sum} / {n}")
    print(f"L_{k}({n}) = {result}")

    print(f"\nTherefore, the number of nonzero coefficients of order {n} is {result}.")

# Set n=10 (order) and k=2 (generators X and Y)
ORDER = 10
GENERATORS = 2
calculate_bch_nonzero_coeffs(ORDER, GENERATORS)