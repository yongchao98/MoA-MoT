import random
from itertools import product

def is_probably_prime(n, k=5):
    """
    Miller-Rabin primality test.
    A probabilistic test that is very accurate for large numbers.
    k is the number of rounds of testing.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    # Write n-1 as 2^r * d
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1

    # Witness loop for k rounds
    for _ in range(k):
        a = random.randrange(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False  # n is composite
    return True  # n is probably prime

def solve_gol_primes_problem():
    """
    This function executes the plan to find the solution.
    1. Generates the first 1,000 primes containing only digits 0 and 1.
    2. Counts how many of these primes, when treated as GoL patterns, would die out.
    """
    # Step 1: Generate the first 1,000 primes (Set C).
    # We generate numbers by their number of digits to ensure we find the first 1,000 in numerical order.
    primes_C = []
    limit = 1000
    num_digits = 2
    
    while len(primes_C) < limit:
        # Generate all binary-like strings of a given length, always starting with '1'.
        for tail in product('01', repeat=num_digits - 1):
            s_num = '1' + ''.join(tail)
            num = int(s_num)
            
            if is_probably_prime(num):
                primes_C.append(num)
                if len(primes_C) == limit:
                    break
        num_digits += 1

    # Step 2: Count how many patterns will die out.
    # This happens if the prime's string representation does not contain "111".
    dying_patterns_count = 0
    for prime in primes_C:
        if "111" not in str(prime):
            dying_patterns_count += 1
            
    # As requested, output the final number.
    print("The number of members of C that will ultimately die out is:")
    print(dying_patterns_count)

solve_gol_primes_problem()