import math

def get_sum_of_distinct_prime_factors(n):
    """
    Calculates the sum of the distinct prime factors of a number.
    For n=1, the sum is 0.
    """
    if n == 1:
        return 0
    
    factors = set()
    d = 2
    temp = n
    
    # Find all factors of 2
    if temp % 2 == 0:
        factors.add(2)
        while temp % 2 == 0:
            temp //= 2
            
    # Find all odd factors
    d = 3
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 2
        
    # If temp is still greater than 1, it must be a prime factor itself
    if temp > 1:
        factors.add(temp)
        
    return sum(factors)

def find_smallest_number():
    """
    Finds the smallest integer N = a + b, where a != b, and the sum of the 
    distinct prime factors of a and b is 20.
    """
    target_sum_of_factors = 20
    n = 1
    
    while True:
        # We check potential sums N starting from 1
        # For each N, we check pairs (a, b) such that a + b = N and a < b.
        # This ensures a and b are distinct and we don't check pairs twice.
        for a in range(1, (n // 2) + 1):
            b = n - a
            
            # Since the loop condition is a <= n/2, a can equal b only if n is even
            # and a = n/2. So, we explicitly skip this case.
            if a == b:
                continue

            # Calculate the sum of prime factors for both numbers
            sopf_a = get_sum_of_distinct_prime_factors(a)
            sopf_b = get_sum_of_distinct_prime_factors(b)
            
            # Check if the sum meets the condition from proposition ¬Q
            if sopf_a + sopf_b == target_sum_of_factors:
                print(f"The number N that satisfies the condition is {n}.")
                print(f"This is because N can be expressed as the sum of a = {a} and b = {b}.")
                print(f"The numbers a and b are different, satisfying condition ¬R.")
                print(f"The sum of the distinct prime factors of a({a}) is {sopf_a}.")
                print(f"The sum of the distinct prime factors of b({b}) is {sopf_b}.")
                print(f"Their total sum is {sopf_a} + {sopf_b} = {target_sum_of_factors}, satisfying condition ¬Q.")
                print(f"Since we searched for N starting from 1, N={n} is the smallest such integer, satisfying condition ¬P.")
                print("\nThe final equation is:")
                print(f"{a} + {b} = {n}")
                return n
        
        n += 1

# Execute the function to find and print the result.
find_smallest_number()