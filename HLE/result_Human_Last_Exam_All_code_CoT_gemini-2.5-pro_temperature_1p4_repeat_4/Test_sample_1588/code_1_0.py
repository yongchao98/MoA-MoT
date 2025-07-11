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

def moebius(n):
    """Calculates the Moebius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0  # n has a squared prime factor
    return (-1)**len(factors) # (-1)^k where k is the number of distinct prime factors

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def calculate_bch_coefficients(n, k):
    """
    Calculates the number of nonzero coefficients of order n for k generators
    using Witt's formula and prints the calculation step-by-step.
    """
    print(f"The number of nonzero BCH coefficients of order n={n} for k={k} generators is given by Witt's formula:")
    print(f"  (1/{n}) * Σ_{{d|{n}}} μ(d) * {k}^({n}/d)")
    print("-" * 30)

    divisors = get_divisors(n)
    print(f"The divisors 'd' of n={n} are: {divisors}\n")
    
    total_sum = 0
    sum_parts = []
    
    print("Calculating each term of the summation Σ μ(d) * k^(n/d):")
    for d in divisors:
        mu_d = moebius(d)
        power = n // d
        k_pow = k**power
        term = mu_d * k_pow
        total_sum += term
        
        print(f"For d = {d}:")
        print(f"  μ({d}) = {mu_d}")
        print(f"  {k}^({n}/{d}) = {k}^{power} = {k_pow}")
        print(f"  Term: μ({d}) * {k_pow} = {term}")
        sum_parts.append(str(term))

    print("-" * 30)
    print("The full summation is:")
    print(f"  Σ μ(d) * k^(n/d) = {' + '.join(sum_parts).replace('+ -', '- ')}")
    print(f"  Sum = {total_sum}")
    print("-" * 30)

    result = total_sum // n
    
    print("Finally, the result is:")
    print(f"  (1/{n}) * {total_sum} = {result}")
    
    print(f"\nThus, the number of nonzero coefficients of order {n} is {result}.")

# Set parameters for the specific problem
order_n = 10
generators_k = 2

# Run the calculation
calculate_bch_coefficients(order_n, generators_k)

print("\n<<<99>>>")