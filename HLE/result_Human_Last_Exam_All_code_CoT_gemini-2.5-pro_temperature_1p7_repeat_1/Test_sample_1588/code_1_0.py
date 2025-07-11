import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors of num."""
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """Computes the Möbius function μ(n)."""
    if n == 1:
        return 1
    
    prime_factors = get_prime_factorization(n)
    
    # Return 0 if there is any squared prime factor
    for p in prime_factors:
        if prime_factors[p] > 1:
            return 0
    
    # Return (-1)^k where k is the number of distinct prime factors
    if len(prime_factors) % 2 != 0:
        return -1
    else:
        return 1

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_coefficients():
    """
    Calculates the number of nonzero coefficients of a given order in the
    BCH expansion using Witt's formula and prints the calculation steps.
    """
    n = 10  # Order of the coefficients
    k = 2   # Number of generators (X, Y)

    divisors = get_divisors(n)
    
    total_sum = 0
    term_strings = []
    value_strings = []

    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term_value = mu_d * (k ** power)
        total_sum += term_value
        
        term_strings.append(f"μ({d})*{k}^({n}/{d})")
        value_strings.append(f"({mu_d})*({k**power})")

    result = total_sum // n
    
    # Print the detailed calculation
    print("The number of nonzero coefficients is given by Witt's formula:")
    print(f"d(n) = (1/n) * Σ_{{d|n}} μ(d) * k^(n/d)\n")
    print(f"For n={n} and k={k}:")
    
    # Build and print the full equation with symbols
    full_equation_sym = f"d({n}) = (1/{n}) * [ " + " + ".join(term_strings) + " ]"
    print(full_equation_sym)
    
    # Build and print the equation with calculated values
    # Replace "+ (-" with "- " for cleaner output
    full_equation_val = f"d({n}) = (1/{n}) * [ " + " + ".join(value_strings) + " ]"
    full_equation_val = full_equation_val.replace("+ (-", "- ")
    print(full_equation_val)

    print(f"d({n}) = (1/{n}) * [ {total_sum} ]")
    print(f"d({n}) = {result}")

# Run the solver
solve_bch_coefficients()
print("\n<<<99>>>")