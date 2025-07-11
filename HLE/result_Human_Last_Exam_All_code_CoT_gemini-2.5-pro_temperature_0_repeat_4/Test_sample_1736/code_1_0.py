import math

def get_sum_of_distinct_prime_factors(n):
    """
    Calculates the sum of distinct prime factors of a given integer n.
    For example, for n=12 (2*2*3), the distinct prime factors are 2 and 3, so the sum is 2+3=5.
    """
    if n < 2:
        return 0
    
    factors = set()
    d = 2
    temp_n = n
    
    # Find all factors of 2
    if temp_n % 2 == 0:
        factors.add(2)
        while temp_n % 2 == 0:
            temp_n //= 2
            
    # Find all odd factors
    d = 3
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 2
        
    # If temp_n is still greater than 1, it must be a prime factor itself
    if temp_n > 1:
        factors.add(temp_n)
        
    return sum(factors)

def find_smallest_N():
    """
    Finds the smallest integer N that is a sum of two different integers a and b,
    such that the sum of their distinct prime factors is 20.
    """
    # We search for N by checking sums in increasing order.
    # The smallest sum of two different positive integers is 1+2=3.
    # We set a reasonable upper limit for the search.
    for N in range(3, 100):
        # For each sum N, check all pairs (a, b) where a+b=N and a < b.
        for a in range(1, N // 2 + 1):
            b = N - a
            
            # Condition: a and b must be different.
            # Our loop a < b (since a <= N/2) ensures this.
            if a == b:
                continue

            # Calculate the sum of prime factors for a and b.
            spf_a = get_sum_of_distinct_prime_factors(a)
            spf_b = get_sum_of_distinct_prime_factors(b)

            # Condition: The sum of prime factors must be 20.
            if spf_a + spf_b == 20:
                # Since we are iterating N in increasing order, the first N we find
                # is guaranteed to be the smallest.
                print(f"Found the smallest N = {N}")
                print(f"The two numbers are a = {a} and b = {b}.")
                print(f"The sum of distinct prime factors of {a} is {spf_a}.")
                print(f"The sum of distinct prime factors of {b} is {spf_b}.")
                print(f"The total sum of prime factors is {spf_a} + {spf_b} = 20.")
                print("\nThe final equation is:")
                print(f"{a} + {b} = {N}")
                return N
    return None

# Execute the search and find the number.
final_N = find_smallest_N()

# The final answer is the number N itself.
# print(f"\n<<<The number that satisfies the condition is {final_N}>>>")
# The above line is commented out to only output the final answer in the required format.
print(f"\n<<<{final_N}>>>")