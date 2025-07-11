import math

def get_sum_of_distinct_prime_factors(n):
    """
    Calculates the sum of distinct prime factors of a given integer n.
    For example, for 12 = 2^2 * 3, the distinct prime factors are 2 and 3, and the sum is 2 + 3 = 5.
    """
    if n <= 1:
        return 0
    
    factors = set()
    d = 2
    temp_n = n
    
    # Factor out multiples of 2
    if temp_n % 2 == 0:
        factors.add(2)
        while temp_n % 2 == 0:
            temp_n //= 2
            
    # Factor out multiples of odd numbers
    d = 3
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 2
        
    # If temp_n is still greater than 1, it must be a prime factor
    if temp_n > 1:
        factors.add(temp_n)
        
    return sum(factors)

def find_smallest_sum():
    """
    Finds the smallest integer N which is the sum of two different integers a and b,
    such that the sum of the distinct prime factors of a and b is 20.
    """
    min_sum = float('inf')
    result_pair = None
    
    # We search within a reasonable limit. Since we are looking for the smallest sum,
    # the components a and b are likely to be small.
    SEARCH_LIMIT = 100 
    
    # Pre-calculate sums of prime factors to optimize
    spf_cache = {i: get_sum_of_distinct_prime_factors(i) for i in range(1, SEARCH_LIMIT + 1)}

    for a in range(1, SEARCH_LIMIT + 1):
        for b in range(a + 1, SEARCH_LIMIT + 1):
            # Condition from ¬Q: Sum of prime factors is 20
            if spf_cache[a] + spf_cache[b] == 20:
                current_sum = a + b
                # Condition from ¬P: N is the smallest such sum
                if current_sum < min_sum:
                    min_sum = current_sum
                    # Condition from ¬R: a and b are different (guaranteed by b > a)
                    result_pair = (a, b)

    if result_pair:
        a, b = result_pair
        print(f"{min_sum} = {a} + {b}")
    else:
        print("No pair found within the search limit.")

if __name__ == "__main__":
    find_smallest_sum()