import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors of num and their powers."""
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

def mobius_mu(n):
    """
    Calculates the Möbius function mu(n).
    - mu(1) = 1
    - mu(n) = 0 if n has a squared prime factor.
    - mu(n) = (-1)^k if n is a product of k distinct primes.
    """
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
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_witt():
    """
    Calculates the number of nonzero coefficients for order 10 in the
    BCH expansion using Witt's formula and prints the full calculation.
    """
    n = 10  # Order of the term
    k = 2   # Number of generators (X and Y)

    # Witt's formula: L_n(k) = (1/n) * sum_{d|n} mu(d) * k^(n/d)
    
    divisors_of_n = get_divisors(n)
    
    sum_terms = 0
    
    # Store parts of the equation string
    equation_parts = []
    
    for d in divisors_of_n:
        mu_d = mobius_mu(d)
        exponent = n // d
        term = mu_d * (k ** exponent)
        sum_terms += term
        
        # Build the string for each term in the sum, e.g., "((-1) * 2^5)"
        equation_parts.append(f"(({mu_d}) * {k}^{exponent})")

    result = sum_terms // n
    
    # Construct the final equation string
    equation_string = f"(1/{n}) * ( " + " + ".join(equation_parts) + f" ) = {result}"
    
    print("The number of nonzero coefficients is calculated using Witt's formula L_n(k) for n=10, k=2.")
    print("The formula is: L_n(k) = (1/n) * sum_{d|n} mu(d) * k^(n/d)")
    print(f"\nThe calculation for n={n}, k={k} and divisors {divisors_of_n} is:")
    print(equation_string)

solve_bch_witt()