import math

def get_prime_factors(num):
    """
    Returns a dictionary of prime factors and their multiplicities.
    e.g., get_prime_factors(12) -> {2: 2, 3: 1}
    """
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
    """Calculates the Moebius function μ(n)."""
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    if n == 1:
        return 1
    
    factors = get_prime_factors(n)
    # If any prime factor has a power greater than 1, μ(n) is 0
    if any(p > 1 for p in factors.values()):
        return 0
    
    # Otherwise, μ(n) is (-1)^k where k is the number of distinct prime factors
    return (-1) ** len(factors)

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_nonzero_coeffs(n, k):
    """
    Calculates and prints the number of nonzero coefficients of order n in the BCH expansion
    using Witt's formula for a free Lie algebra with k generators.
    """
    divs = get_divisors(n)
    
    total_sum = 0
    symbolic_parts = []
    calculated_parts = []
    
    for d in divs:
        mu_d = moebius(d)
        power = n // d
        base_power_val = k ** power
        term_val = mu_d * base_power_val
        total_sum += term_val
        
        # Build symbolic part for the equation, e.g., "(-1)*(2^5)"
        symbolic_parts.append(f"({mu_d})*({k}^{power})")
        
        # Build calculated part, e.g., "- 32"
        if term_val >= 0:
            calculated_parts.append(f"+ {term_val}")
        else:
            calculated_parts.append(f"- {-term_val}")

    result = total_sum // n
    
    # Adjust the first calculated part to remove the leading "+"
    if calculated_parts[0].startswith('+ '):
        calculated_parts[0] = calculated_parts[0][2:]

    # Construct the final equation string for printing
    symbolic_str = " + ".join(symbolic_parts)
    calculated_str = " ".join(calculated_parts)
    
    print(f"The number of nonzero coefficients of order {n} is: {result}")
    print("\nThis is found using Witt's formula L_n(k) = (1/n) * Σ_{d|n} μ(d) * k^(n/d)")
    print(f"For n={n} and k={2}, the divisors are {divs}.")
    print("\nThe full calculation is:")
    print(f"(1/{n}) * ( {symbolic_str} )")
    print(f"= (1/{n}) * ( {calculated_str} )")
    print(f"= (1/{n}) * ({total_sum})")
    print(f"= {result}")

if __name__ == '__main__':
    # Order of the expansion
    order_n = 10
    # Number of generators (X and Y)
    generators_k = 2
    solve_bch_nonzero_coeffs(order_n, generators_k)
<<<99>>>