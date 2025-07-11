import math

def sum_prime_factors(n):
    """
    Calculates the sum of the unique prime factors of a given number n.
    For example, sum_prime_factors(12) returns 2 + 3 = 5.
    """
    if n <= 1:
        return 0
    
    factors = set()
    d = 2
    temp_n = n
    
    # Find factors up to the square root of n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
        
    # If there's a remainder, it must be a prime factor
    if temp_n > 1:
        factors.add(temp_n)
        
    return sum(factors)

def find_the_number():
    """
    Finds the smallest integer N that is the sum of two different integers a and b,
    where the sum of the prime factors of a and b is 20.
    """
    
    # We are looking for the smallest N = a + b satisfying the conditions.
    # We can iterate through N starting from the smallest possible sum of two
    # different positive integers (1 + 2 = 3).
    N = 3
    while True:
        # For a given N, check all pairs (a, b) that sum to N.
        # We only need to check a up to N/2 to avoid duplicate pairs.
        for a in range(1, (N // 2) + 1):
            b = N - a

            # Condition ¬R: a and b must be different.
            # This loop structure ensures a <= b, so we just need to skip the a == b case.
            if a == b:
                continue

            # Condition ¬Q: The sum of the prime factors of a and b must be 20.
            spf_a = sum_prime_factors(a)
            spf_b = sum_prime_factors(b)
            
            if spf_a + spf_b == 20:
                # Condition ¬P: We need the smallest N. Since we are iterating N
                # upwards, the first N that satisfies the conditions is the smallest.
                print("Found the smallest number N that satisfies the conditions.")
                print(f"The two different integers are a = {a} and b = {b}.")
                print(f"The sum of prime factors of {a} is {spf_a}.")
                print(f"The sum of prime factors of {b} is {spf_b}.")
                print(f"Their total sum of prime factors is {spf_a + spf_b}.")
                print("\nThe final equation is:")
                print(f"{a} + {b} = {N}")
                return N
        
        N += 1

if __name__ == "__main__":
    find_the_number()
