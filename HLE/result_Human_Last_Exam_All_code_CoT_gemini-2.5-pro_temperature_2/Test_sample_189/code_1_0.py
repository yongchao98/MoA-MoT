import math

def find_prime_chain():
    """
    Finds the largest prime p < 1,000,000 of the form p=4u+1, where u is
    a prime of the form u=4v+1, and v is also a prime.
    """
    limit = 1_000_000

    # Step 1: Use Sieve of Eratosthenes to find all primes up to the limit.
    # is_prime[i] will be True if i is prime, and False otherwise.
    is_prime = [True] * limit
    is_prime[0] = is_prime[1] = False
    for start_num in range(2, int(math.sqrt(limit)) + 1):
        if is_prime[start_num]:
            # Mark all multiples of start_num as not prime
            for multiple in range(start_num * start_num, limit, start_num):
                is_prime[multiple] = False

    # Step 2: Iterate downwards from the limit to find the largest p.
    # The first valid chain found will be the one with the largest p.
    for p in range(limit - 1, 1, -1):
        # Condition 1: p must be prime
        if is_prime[p]:
            # Condition 2: p must be of the form 4u+1
            if (p - 1) % 4 == 0:
                u = (p - 1) // 4
                
                # Condition 3: u must be prime
                if u > 0 and is_prime[u]:
                    # Condition 4: u must be of the form 4v+1
                    if (u - 1) % 4 == 0:
                        v = (u - 1) // 4
                        
                        # Condition 5: v must be prime
                        if v > 0 and is_prime[v]:
                            # Success! We found the chain.
                            print("Found the largest prime p in the chain.")
                            print(f"The equation for p is: {p} = 4 * {u} + 1")
                            print(f"The equation for u is: {u} = 4 * {v} + 1")
                            print(f"The full chain (p:u:v) is: {p}:{u}:{v}")
                            return

    print("No such prime chain was found within the limit.")

if __name__ == '__main__':
    find_prime_chain()