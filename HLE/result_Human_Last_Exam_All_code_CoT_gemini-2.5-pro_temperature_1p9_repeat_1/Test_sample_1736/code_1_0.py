import math

def get_sum_of_prime_factors(n, cache):
    """
    Calculates the sum of the unique prime factors of a number n.
    Uses a cache to store results for previously computed numbers.
    e.g., for 12 = 2*2*3, unique prime factors are {2, 3}, sum is 5.
    """
    if n in cache:
        return cache[n]
    
    if n < 2:
        return 0
    
    original_n = n
    factors = set()
    
    # Handle the factor of 2
    while n % 2 == 0:
        factors.add(2)
        n //= 2
    
    # Handle odd factors up to the square root of n
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.add(i)
            n //= i
            
    # If n remains as a number greater than 2, it is a prime factor
    if n > 2:
        factors.add(n)
        
    result = sum(factors)
    cache[original_n] = result
    return result

def solve_number_puzzle():
    """
    Solves the logical puzzle by finding the smallest N.
    """
    min_N = float('inf')
    best_a = None
    best_b = None
    
    # A cache to memoize results for the helper function
    prime_factor_sum_cache = {}

    # A search limit of 100 is sufficient, as the numbers involved are small.
    SEARCH_LIMIT = 100

    print("We need to find the smallest N = a + b where:")
    print("1. a and b are different numbers.")
    print("2. The sum of the unique prime factors of a and b is 20.")
    print("\nSearching for the pair (a, b) that satisfies these conditions and minimizes their sum N...")

    for a in range(2, SEARCH_LIMIT):
        sum_pf_a = get_sum_of_prime_factors(a, prime_factor_sum_cache)
        
        # Start b from a + 1 to ensure a != b and avoid duplicate pair checks
        for b in range(a + 1, SEARCH_LIMIT):
            sum_pf_b = get_sum_of_prime_factors(b, prime_factor_sum_cache)
            
            # Check if the sum of prime factors is 20
            if sum_pf_a + sum_pf_b == 20:
                current_N = a + b
                
                # Check if this sum N is the smallest one found so far
                if current_N < min_N:
                    min_N = current_N
                    best_a = a
                    best_b = b

    # Print the details of the solution found
    if best_a is not None:
        sum_a = get_sum_of_prime_factors(best_a, prime_factor_sum_cache)
        sum_b = get_sum_of_prime_factors(best_b, prime_factor_sum_cache)
        
        print(f"\nOptimal pair found: a = {best_a}, b = {best_b}.")
        print(f"The sum of prime factors of {best_a} is {sum_a}.")
        print(f"The sum of prime factors of {best_b} is {sum_b}.")
        print(f"Total sum of prime factors = {sum_a} + {sum_b} = {sum_a + sum_b}.")

        print("\nThe smallest number N satisfying the conditions is the sum of this pair.")
        print(f"Final equation: {best_a} + {best_b} = {min_N}")
    else:
        print("No solution found within the search limit.")

solve_number_puzzle()