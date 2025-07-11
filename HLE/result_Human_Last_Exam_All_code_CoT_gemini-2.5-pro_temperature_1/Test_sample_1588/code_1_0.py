import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of num.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1}
    """
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
    """
    Calculates the Mobius function mu(n).
    - mu(1) = 1
    - mu(n) = 0 if n has a squared prime factor.
    - mu(n) = (-1)^k if n is a product of k distinct primes.
    """
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0 # Has a squared prime factor
            
    # n is square-free, return (-1)^k where k is the number of distinct prime factors
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a sorted list of positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_coeffs(n, k):
    """
    Calculates the number of nonzero coefficients of order n in BCH expansion
    for k generators using Witt's formula.
    """
    divisors = get_divisors(n)
    
    total_sum = 0
    
    print(f"Calculating the number of nonzero coefficients of order n={n} for k={k} generators.")
    print("Using Witt's Formula: L_k(n) = (1/n) * sum_{d|n} mu(d) * k^(n/d)")
    
    # Build the formula string for display
    formula_parts = []
    value_parts = []
    
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term = mu_d * (k ** power)
        total_sum += term
        
        formula_parts.append(f"μ({d})*{k}^({n}/{d})")
        value_parts.append(f"({mu_d})*({k**power})")
        
    formula_str = " + ".join(formula_parts)
    values_str = " + ".join(value_parts)
    
    print(f"\nL_{k}({n}) = (1/{n}) * ( {formula_str} )")
    print(f"       = (1/{n}) * ( {values_str} )")
    print(f"       = (1/{n}) * ( {total_sum} )")
    
    result = total_sum // n
    
    print(f"       = {result}")
    print(f"\nThe number of nonzero coefficients of order {n} is {result}.")

# Main execution
if __name__ == "__main__":
    order_n = 10
    generators_k = 2
    solve_bch_coeffs(order_n, generators_k)