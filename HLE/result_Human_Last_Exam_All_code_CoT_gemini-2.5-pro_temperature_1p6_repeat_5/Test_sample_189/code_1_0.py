import math

def find_largest_chained_prime():
    """
    Finds the largest prime p < 1,000,000 such that p = 4u + 1 and u = 4v + 1,
    where p, u, and v are all prime numbers.
    The architecture of the ancient computer supports integers up to 999,999.
    """
    limit = 999999

    # Step 1: Use Sieve of Eratosthenes to find all primes up to the limit.
    # is_prime[i] will be true if i is prime.
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.sqrt(limit)) + 1):
        if is_prime[i]:
            for multiple in range(i*i, limit + 1, i):
                is_prime[multiple] = False

    # Step 2: Determine the search range for v.
    # p = 4u + 1 = 4(4v + 1) + 1 = 16v + 5
    # To find the largest p, we start with the largest possible v.
    # 16v + 5 <= limit  => v <= (limit - 5) / 16
    max_v = (limit - 5) // 16

    # Step 3: Iterate v downwards to find the largest p first.
    for v in range(max_v, 1, -1):
        # Check if v is prime
        if is_prime[v]:
            u = 4 * v + 1
            # u must be within the bounds of our sieve
            if u > limit:
                continue

            # Check if u is prime
            if is_prime[u]:
                p = 4 * u + 1
                # p must be within the bounds of our sieve
                if p > limit:
                    continue
                
                # Check if p is prime
                if is_prime[p]:
                    # Since we are iterating downwards, the first result found is the largest.
                    # Output each number in the final chain as requested.
                    print(f"Largest prime p found: {p}")
                    print(f"The equation is: {p} = 4 * {u} + 1")
                    print(f"Where u = {u}, which is derived from: u = 4 * {v} + 1")
                    print(f"The final chain in p:u:v format is:")
                    print(f"{p}:{u}:{v}")
                    return

    print("No such prime chain found within the given limit.")

if __name__ == "__main__":
    find_largest_chained_prime()