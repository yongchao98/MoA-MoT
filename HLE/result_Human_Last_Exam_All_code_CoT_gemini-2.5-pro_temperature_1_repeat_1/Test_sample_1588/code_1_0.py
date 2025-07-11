import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of a number and their powers.
    Example: get_prime_factorization(12) -> {2: 2, 3: 1}
    """
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

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    mu(n) = 1 if n is a square-free positive integer with an even number of prime factors.
    mu(n) = -1 if n is a square-free positive integer with an odd number of prime factors.
    mu(n) = 0 if n has a squared prime factor.
    """
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
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def main():
    """
    Calculates the number of nonzero coefficients of a given order in the
    BCH expansion using Witt's formula and prints the calculation steps.
    """
    n = 10  # The order
    k = 2   # The number of generators (X and Y)

    divisors = get_divisors(n)
    
    terms = []
    for d in divisors:
        mu_d = mobius(d)
        power_val = k**(n // d)
        term_val = mu_d * power_val
        terms.append(term_val)
        
    total_sum = sum(terms)
    result = total_sum // n
    
    # Build the string for the sum part of the equation for printing
    sum_str_parts = []
    for term in terms:
        if not sum_str_parts: # First term
            sum_str_parts.append(str(term))
        elif term >= 0:
            sum_str_parts.append(f"+ {term}")
        else:
            sum_str_parts.append(f"- {abs(term)}")
    
    sum_str = " ".join(sum_str_parts)

    print(f"The number of nonzero coefficients is L_k(n) for n={n}, k={k}.")
    print("The calculation using Witt's formula is:")
    print(f"L_2(10) = (1/{n}) * ({sum_str})")
    print(f"       = (1/{n}) * {total_sum}")
    print(f"       = {result}")

if __name__ == "__main__":
    main()