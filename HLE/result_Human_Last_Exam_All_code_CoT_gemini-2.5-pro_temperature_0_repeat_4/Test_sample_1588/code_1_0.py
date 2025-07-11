import math

def get_prime_factorization(num):
    """Returns the prime factorization of a number as a dictionary."""
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
    """Calculates the Möbius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0  # n has a squared prime factor, so mu(n) is 0
    # If square-free, mu(n) is (-1)^k where k is the number of distinct prime factors
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_coefficients():
    """
    Calculates the number of nonzero BCH coefficients for a given order n
    and number of generators k, and prints the detailed calculation.
    """
    n = 10
    k = 2

    divisors = get_divisors(n)
    total_sum = 0

    # Build the string for the first line of the equation
    # e.g., (1/10) * (1 * 2^10 - 1 * 2^5 - ...)
    sum_str_parts = []
    for d in divisors:
        mu_d = mobius(d)
        if mu_d == 0:
            continue # Skip terms where mu(d) is 0
        power = n // d
        
        sign = " + " if mu_d > 0 else " - "
        if not sum_str_parts: # First term logic
            sign = "" if mu_d > 0 else "- "
            
        sum_str_parts.append(f"{sign}{abs(mu_d)} * {k}^{power}")

    # Build the string for the second line of the equation (with evaluated powers)
    # e.g., (1/10) * (1024 - 32 - 4 + 2)
    sum_val_parts = []
    for d in divisors:
        mu_d = mobius(d)
        if mu_d == 0:
            continue
        power = n // d
        term_val = k**power
        total_sum += mu_d * term_val
        
        sign = " + " if mu_d > 0 else " - "
        if not sum_val_parts: # First term logic
            sign = "" if mu_d > 0 else "- "
            
        sum_val_parts.append(f"{sign}{term_val}")

    final_result = total_sum // n

    # Print the step-by-step calculation
    print(f"The number of nonzero coefficients of order {n} is given by Witt's formula L_n(k).")
    print(f"For n={n} and k={k}, the divisors of {n} are {divisors}.")
    print(f"The corresponding values of the Möbius function mu(d) are {[mobius(d) for d in divisors]}.\n")
    
    print("The calculation is as follows:")
    # Line 1: Symbolic form with values
    print(f"L_{n}({k}) = (1/{n}) * ({''.join(sum_str_parts)})")
    # Line 2: Evaluated powers
    print(f"       = (1/{n}) * ({''.join(sum_val_parts)})")
    # Line 3: Sum evaluated
    print(f"       = (1/{n}) * ({total_sum})")
    # Line 4: Final result
    print(f"       = {final_result}")

if __name__ == "__main__":
    solve_bch_coefficients()