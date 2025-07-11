import math

def get_sum_of_prime_factors(n):
    """
    Calculates the sum of the unique prime factors of a given number n.
    For example, for 12 (2*2*3), the unique prime factors are {2, 3}, and the sum is 5.
    """
    if n < 2:
        return 0
    factors = set()
    temp_n = n
    
    # Handle the factor 2
    if temp_n % 2 == 0:
        factors.add(2)
        while temp_n % 2 == 0:
            temp_n //= 2
            
    # Handle odd factors
    for i in range(3, int(math.sqrt(temp_n)) + 1, 2):
        if temp_n % i == 0:
            factors.add(i)
            while temp_n % i == 0:
                temp_n //= i
                
    # If temp_n is a prime number greater than 2 after division
    if temp_n > 2:
        factors.add(temp_n)
        
    return sum(factors)

def solve_puzzle():
    """
    Solves the logic puzzle by finding the smallest integer N that satisfies the conditions.

    The negated propositions require us to find the smallest N such that:
    1. N = a + b, where a and b are different integers (from ¬P and ¬R).
    2. The sum of the prime factors of a and b is 20 (from ¬Q).
    """
    
    # We will search for N by checking sums in increasing order.
    # The smallest possible N from two distinct integers >= 2 is 2+3=5.
    n = 4
    spf_cache = {}

    while True:
        n += 1
        # To find pairs (a, b) with a+b=n and a < b, we can iterate a from 2 up to n/2.
        for a in range(2, (n // 2) + 1):
            if a * 2 == n:  # This ensures a and b are different
                continue

            b = n - a

            # Use a cache for performance
            if a not in spf_cache:
                spf_cache[a] = get_sum_of_prime_factors(a)
            if b not in spf_cache:
                spf_cache[b] = get_sum_of_prime_factors(b)
            
            spf_a = spf_cache[a]
            spf_b = spf_cache[b]

            if spf_a + spf_b == 20:
                # We found the smallest N because we are iterating N in increasing order.
                print(f"{n} = {a} + {b}")
                return n

# Execute the solver to find and print the answer.
found_n = solve_puzzle()
# The required output format is just the answer content.
# The code above prints the equation, so this final print wraps the number.
print(f"<<<{found_n}>>>")
