import math

def get_prime_factorization(num):
    """Returns the prime factorization of num as a dictionary."""
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

def moebius(n):
    """Computes the Moebius function μ(n)."""
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
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

def solve_bch_coeffs(n, k):
    """
    Calculates the number of nonzero coefficients of order n in the BCH expansion
    for k generators using Witt's formula.
    """
    print(f"Calculating the number of nonzero coefficients for order n={n} and k={k} generators.")
    print("Using Witt's formula: d(n) = (1/n) * Σ_{d|n} μ(n/d) * k^d\n")

    divs = get_divisors(n)
    print(f"The divisors 'd' of n={n} are: {divs}\n")

    total_sum = 0
    sum_terms = []

    print("Calculating each term in the summation Σ μ(n/d) * k^d:")
    for d in divs:
        n_over_d = n // d
        mu_val = moebius(n_over_d)
        k_power_d = k**d
        term = mu_val * k_power_d
        
        print(f"For d = {d}:")
        print(f"  μ(n/d) * k^d = μ({n}/{d}) * {k}^{d}")
        print(f"               = μ({n_over_d}) * {k_power_d}")
        print(f"               = {mu_val} * {k_power_d}")
        print(f"               = {term}")

        total_sum += term
        sum_terms.append(str(term))

    print("\nThe full summation is:")
    sum_expression = " + ".join(sum_terms).replace("+ -", "- ")
    print(f"Σ μ(n/d) * k^d = {sum_expression} = {total_sum}\n")
    
    if total_sum % n != 0:
         print(f"Warning: Sum ({total_sum}) is not divisible by n ({n}). Check logic.")
         result = total_sum / n
    else:
        result = total_sum // n

    print(f"The final result is d({n}) = (1/{n}) * {total_sum} = {result}")

# Parameters for the problem
order_n = 10
generators_k = 2

solve_bch_coeffs(order_n, generators_k)