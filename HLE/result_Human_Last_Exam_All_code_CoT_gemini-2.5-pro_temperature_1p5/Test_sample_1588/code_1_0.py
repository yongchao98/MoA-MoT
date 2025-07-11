import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of num.
    Returns None if a factor is repeated (e.g., for 12 = 2*2*3).
    """
    factors = {}
    d = 2
    temp_num = num
    while d * d <= temp_num:
        if temp_num % d == 0:
            if d in factors:
                return None # Repeated factor
            factors[d] = 1
            temp_num //= d
            # Check if the factor is repeated
            if temp_num % d == 0:
                return None
        d += 1
    if temp_num > 1:
        if temp_num in factors:
            return None # Repeated factor
        factors[temp_num] = 1
    return factors

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    mu(n) = 1 if n is a square-free positive integer with an even number of prime factors.
    mu(n) = -1 if n is a square-free positive integer with an odd number of prime factors.
    mu(n) = 0 if n has a squared prime factor.
    mu(n) = 1 for n=1.
    """
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    
    if factors is None: # Has a squared prime factor
        return 0
    
    # Square-free, return (-1)^k where k is the number of prime factors
    return (-1) ** len(factors)

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_coeffs():
    """
    Calculates the number of nonzero coefficients for a given order in the BCH expansion.
    """
    n = 10  # Order of the expansion
    k = 2   # Number of non-commuting variables (X and Y)

    divisors = get_divisors(n)
    
    total_sum = 0
    term_strings = []
    value_strings = []
    
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term = mu_d * (k ** power)
        total_sum += term
        
        term_strings.append(f"mu({d})*({k}^{power})")
        value_strings.append(f"{mu_d}*({k**power})")

    result = total_sum // n

    print(f"The number of nonzero coefficients is given by Witt's formula L_n(k) = (1/n) * sum_{{d|n}} mu(d) * k^(n/d)")
    print(f"For n={n} and k={k}, the divisors of {n} are: {divisors}")
    print("\nCalculation steps:")
    print(f"L_{n}({k}) = (1/{n}) * ( " + " + ".join(term_strings) + " )")
    print(f"L_{n}({k}) = (1/{n}) * ( " + " + ".join(value_strings).replace("+-", "- ") + " )")
    print(f"L_{n}({k}) = (1/{n}) * ( {total_sum} )")
    print(f"L_{n}({k}) = {result}")
    
    print(f"\nThe number of nonzero coefficients of order {n} is {result}.")

if __name__ == "__main__":
    solve_bch_coeffs()
    print("<<<99>>>")
