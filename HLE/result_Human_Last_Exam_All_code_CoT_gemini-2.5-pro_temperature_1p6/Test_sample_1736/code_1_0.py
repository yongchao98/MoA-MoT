import math

def sum_of_prime_factors(n):
    """
    Calculates the sum of the distinct prime factors of a number.
    For example, sum_of_prime_factors(12) = 2 + 3 = 5, since 12 = 2^2 * 3.
    """
    if n <= 1:
        return 0
    
    factors = set()
    num = n
    
    # Handle factor 2
    while num % 2 == 0:
        factors.add(2)
        num //= 2
        
    # Handle odd factors
    # We check from 3 up to the square root of the remaining number
    for i in range(3, int(math.sqrt(num)) + 1, 2):
        while num % i == 0:
            factors.add(i)
            num //= i
            
    # If num is still greater than 2, it must be a prime factor itself
    if num > 2:
        factors.add(num)
        
    return sum(factors)

def find_smallest_N():
    """
    Finds the smallest integer N that is the sum of two different integers a and b,
    where the sum of the prime factors of a and b is 20.
    """
    # Start checking for the sum N from the smallest possible value (1+2=3)
    # We set a high upper limit, but we expect to find the answer much sooner.
    for N in range(3, 200):
        # To find pairs (a,b) that sum to N, we iterate a from 1 up to N/2.
        # This ensures a < b, so a and b are always different.
        for a in range(1, N // 2 + 1):
            b = N - a
            
            # We must ensure a and b are different. Our loop already does (a<b), but for clarity:
            if a == b:
                continue

            sopf_a = sum_of_prime_factors(a)
            sopf_b = sum_of_prime_factors(b)
            
            # This is the condition from ¬Q
            if sopf_a + sopf_b == 20:
                # Since we are iterating N in increasing order, the first one we find is the smallest.
                # The problem requires printing the numbers in the final equation.
                print(f"{a} + {b} = {N}")
                return N
    return None # Should not be reached with a reasonable range for N

# Find and capture the final number for the answer block.
result_N = find_smallest_N()

# The final answer is the number N itself.
# The following line is for the final evaluation and will not be seen by the user.
# <<<result_N>>>