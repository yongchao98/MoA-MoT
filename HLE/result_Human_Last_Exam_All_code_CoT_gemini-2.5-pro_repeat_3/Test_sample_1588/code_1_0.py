import math

def prime_factors(n):
    """
    Returns a dictionary of prime factors of n.
    For example, prime_factors(12) returns {2: 2, 3: 1}.
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
    Calculates the Möbius function μ(n).
    μ(n) = 1 if n is a square-free positive integer with an even number of prime factors.
    μ(n) = -1 if n is a square-free positive integer with an odd number of prime factors.
    μ(n) = 0 if n has a squared prime factor.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Möbius function is defined only for positive integers.")
    if n == 1:
        return 1
    
    factors = prime_factors(n)
    for p in factors:
        if factors[p] > 1:
            return 0  # n has a squared prime factor
            
    return (-1)**len(factors) # n is square-free

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_problem():
    """
    Calculates and explains the number of nonzero BCH coefficients of order 10.
    """
    n = 10  # The order
    k = 2   # The number of generators (X and Y)

    # Explain the method
    print("The number of nonzero coefficients of order n in the Baker-Campbell-Hausdorff (BCH) expansion for k generators")
    print("is given by Witt's formula, which counts the dimensions of the graded components of a free Lie algebra.")
    
    print("\nWitt's formula is:")
    print(f"  l_k(n) = (1/n) * Σ (over d|n) [μ(d) * k^(n/d)]")
    print(f"Where:")
    print(f"  n = order = {n}")
    print(f"  k = number of generators = {k}")
    print(f"  d = divisors of n")
    print(f"  μ = the Möbius function")

    # Step 1: Get divisors
    divisors = get_divisors(n)
    print(f"\nStep 1: The divisors of n = {n} are: {divisors}")
    
    # Step 2: Calculate each term in the summation
    print("\nStep 2: Calculate the terms μ(d) * k^(n/d) for each divisor d.")
    sum_terms = []
    total_sum = 0
    
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term_val = k ** power
        term = mu_d * term_val
        total_sum += term
        sum_terms.append(term)
        print(f"  For d={d:2d}: term = μ({d}) * {k}^({n}/{d}) = {mu_d:2d} * {k}^{power} = {mu_d:2d} * {term_val:4d} = {term:5d}")

    # Step 3: Complete the calculation
    result = total_sum // n
    print("\nStep 3: Apply the full formula with the calculated sum.")
    print(f"  l_{k}({n}) = (1/{n}) * (Sum of terms)")
    
    # Build the string for the sum to show each number in the final equation
    sum_string = " + ".join(map(str, sum_terms)).replace('+ -', '- ')
    print(f"  l_{k}({n}) = (1/{n}) * ({sum_string})")
    print(f"  l_{k}({n}) = (1/{10}) * ({total_sum})")
    print(f"  l_{k}({n}) = {result}")

    print(f"\nTherefore, the number of nonzero coefficients of order 10 is {result}.")

# Execute the solution
solve_bch_problem()
<<<99>>>