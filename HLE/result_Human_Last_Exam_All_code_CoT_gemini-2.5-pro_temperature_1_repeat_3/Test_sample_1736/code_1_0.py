import math

def sum_prime_factors(n):
    """
    Calculates the sum of the prime factors of a number (with multiplicity).
    For example, sum_prime_factors(12) = 2 + 2 + 3 = 7.
    """
    if n <= 1:
        return 0
    
    total_sum = 0
    d = 2
    temp_n = n
    
    # Divide by 2 until it's odd
    while temp_n % 2 == 0:
        total_sum += 2
        temp_n //= 2
        
    # Check for odd factors from 3 upwards
    d = 3
    while d * d <= temp_n:
        while temp_n % d == 0:
            total_sum += d
            temp_n //= d
        d += 2
        
    # If temp_n is still a prime number greater than 1
    if temp_n > 1:
        total_sum += temp_n
        
    return total_sum

def find_smallest_sum():
    """
    Finds the smallest integer N which is the sum of two distinct integers a and b,
    where the sum of the prime factors of a and b is 20.
    """
    min_sum = float('inf')
    best_a = -1
    best_b = -1
    
    # A search limit of 100 is sufficient, as the smallest N is found with a and b well below this.
    limit = 100
    spf_cache = {i: sum_prime_factors(i) for i in range(1, limit + 1)}

    # Iterate through all pairs of distinct numbers (a, b)
    for a in range(2, limit):
        for b in range(a + 1, limit): # b > a ensures a and b are different
            
            # Check if the sum of prime factors is 20
            if spf_cache.get(a, 0) + spf_cache.get(b, 0) == 20:
                current_sum = a + b
                # If we find a new smaller sum, update the minimum
                if current_sum < min_sum:
                    min_sum = current_sum
                    best_a = a
                    best_b = b

    # Print the final equation
    if best_a != -1:
        # We need to present the equation. Let's make sure a is the smaller number.
        final_a = min(best_a, best_b)
        final_b = max(best_a, best_b)
        print(f"The number N is the result of the equation: {final_a} + {final_b} = {min_sum}")
    else:
        print("No solution found within the specified limit.")

# Execute the function to find and print the result.
find_smallest_sum()