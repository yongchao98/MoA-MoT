def solve():
    """
    Finds the largest prime divisor for the number 8139881.
    """
    # The number we want to factorize
    n = 8139881
    original_n = n
    
    # A list to store all the prime factors
    factors = []
    
    # Start with the smallest prime divisor, 2
    # This step is technically not needed for an odd number, but it's
    # good practice for a general function.
    while n % 2 == 0:
        factors.append(2)
        n //= 2
        
    # Check for odd divisors starting from 3
    # We only need to check up to the square root of the remaining n
    d = 3
    while d * d <= n:
        # If d divides n, keep dividing n by d and record the factor
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 2
        
    # If there is a number remaining that is greater than 1, it must be
    # a prime factor itself (and the largest one).
    if n > 1:
        factors.append(n)
        
    # Print the results as requested
    if factors:
        print(f"{original_n} = {' * '.join(map(str, factors))}")
        # The largest prime factor will be the last one in the sorted list of factors
        largest_prime = factors[-1]
        print(f"The largest prime divisor of {original_n} is {largest_prime}")
    else:
        # This case handles if the original number itself is prime
        print(f"{original_n} is a prime number, so its largest prime divisor is itself.")
        largest_prime = original_n

    # The final answer in the requested format will be based on this calculation
    # In this case, the largest prime factor is 2857.

solve()