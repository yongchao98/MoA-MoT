import math

def get_prime_factorization(num):
    """Returns the prime factorization of a number as a dictionary of {prime: exponent}."""
    factors = {}
    d = 2
    temp_num = num
    while d * d <= temp_num:
        while (temp_num % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_num //= d
        d += 1
    if temp_num > 1:
       factors[temp_num] = factors.get(temp_num, 0) + 1
    return factors

def mobius(n):
    """Calculates the Möbius function mu(n)."""
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    
    # If any prime factor is squared, mu(n) is 0
    for p in factors:
        if factors[p] > 1:
            return 0
            
    # Otherwise, it's (-1)^k where k is the number of distinct prime factors
    if len(factors) % 2 == 1:
        return -1
    else:
        return 1

def get_divisors(n):
    """Returns a sorted list of all divisors of a number."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def calculate_bch_coefficients():
    """
    Calculates the number of nonzero coefficients for a given order in the
    BCH expansion using Witt's formula and prints the detailed calculation.
    """
    n = 10  # Order
    k = 2   # Number of generators (X and Y)

    print(f"Calculating the number of nonzero BCH coefficients of order n={n} for k={k} generators.")
    print("Using Witt's formula: d(n) = (1/n) * sum over d|n of (mu(d) * k^(n/d))")
    
    divisors = get_divisors(n)
    print(f"\nThe divisors 'd' of n={n} are: {divisors}\n")
    
    terms = []
    total_sum = 0
    
    # Calculate each term in the summation
    for d in divisors:
        mu_d = mobius(d)
        power_val = k**(n // d)
        term = mu_d * power_val
        terms.append(term)
        total_sum += term
        
    print("The final equation and its step-by-step evaluation is:")
    
    # 1. Print the formula with symbols
    sum_formula_parts = []
    for d in divisors:
        sum_formula_parts.append(f"mu({d})*({k}^({n}/{d}))")
    sum_formula_str = " + ".join(sum_formula_parts)
    print(f"Result = ( {sum_formula_str} ) / {n}")

    # 2. Print the formula with numerical values for mu(d)
    sum_eval_parts = []
    for d in divisors:
        sum_eval_parts.append(f"({mobius(d)})*({k}^{n//d})")
    sum_eval_str = " + ".join(sum_eval_parts)
    print(f"Result = ( {sum_eval_str} ) / {n}")

    # 3. Print the formula with evaluated powers
    sum_power_eval_parts = []
    for d in divisors:
        sum_power_eval_parts.append(f"({mobius(d)})*({k**(n//d)})")
    sum_power_eval_str = " + ".join(sum_power_eval_parts)
    print(f"Result = ( {sum_power_eval_str} ) / {n}")
    
    # 4. Print the formula with evaluated terms
    sum_term_vals_str = " + ".join(map(str, terms)).replace('+ -', '- ')
    print(f"Result = ( {sum_term_vals_str} ) / {n}")
    
    # 5. Print the sum
    print(f"Result = {total_sum} / {n}")
    
    # 6. Print the final result
    final_answer = total_sum // n
    print(f"Result = {final_answer}")

if __name__ == '__main__':
    calculate_bch_coefficients()