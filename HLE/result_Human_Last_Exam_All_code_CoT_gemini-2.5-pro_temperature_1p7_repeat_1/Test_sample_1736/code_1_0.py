import math

def get_sum_of_prime_factors(n):
    """
    Calculates the sum of the unique prime factors of a given integer n.
    For n=1, the sum is 0.
    """
    if n <= 1:
        return 0
    
    factors = set()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    if temp_n > 1:
        factors.add(temp_n)
        
    return sum(factors)

def find_smallest_sum():
    """
    Finds the smallest integer N that is the sum of two different integers a and b,
    where the sum of the prime factors of a and b is 20.
    """
    min_sum = float('inf')
    result_a = -1
    result_b = -1
    
    # A search limit of 100 is sufficient to find the smallest pairs.
    search_limit = 100
    
    # Pre-calculate sums of prime factors for efficiency
    spf_cache = {i: get_sum_of_prime_factors(i) for i in range(1, search_limit)}

    for a in range(1, search_limit):
        # We start b from a + 1 to ensure a and b are different.
        for b in range(a + 1, search_limit):
            
            # Condition from ¬Q: The sum of the prime factors of a and b is equal to 20.
            if spf_cache[a] + spf_cache[b] == 20:
                current_sum = a + b
                
                # Condition from ¬P: N is the smallest such sum.
                if current_sum < min_sum:
                    min_sum = current_sum
                    result_a = a
                    result_b = b

    if result_a != -1:
        print("The condition is met by several pairs, the smallest sum N is produced by a pair such as:")
        # We want to output each number in the final equation.
        print(f"{result_a} + {result_b} = {min_sum}")
        print(f"The number N that satisfies the given conditions is: {min_sum}")
    else:
        print("No pair found within the search limit.")

if __name__ == '__main__':
    find_smallest_sum()