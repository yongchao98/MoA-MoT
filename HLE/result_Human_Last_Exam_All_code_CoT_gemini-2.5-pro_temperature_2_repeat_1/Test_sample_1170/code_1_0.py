import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def get_prime_factors(n):
    """Returns a set of unique prime factors of a number."""
    factors = set()
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
    return factors

def solve():
    """
    Solves the problem based on the numerological hypothesis.
    """
    max_object_index = 20

    # Part 1: Find which object indices are prime numbers.
    prime_generators = [i for i in range(1, max_object_index + 1) if is_prime(i)]
    
    # Part 2: Count numbers <= 20 that are powers of a prime but not prime.
    # A number n is a power of a prime if it has only one unique prime factor.
    non_generated_count = 0
    for i in range(1, max_object_index + 1):
        if not is_prime(i) and i > 1:
            factors = get_prime_factors(i)
            if len(factors) == 1:
                non_generated_count += 1
                
    # Format the final output string
    result_string = ",".join(map(str, prime_generators)) + ";" + str(non_generated_count)
    print(result_string)

solve()