import math

def find_special_prime():
    """
    This script searches for a prime number n and an integer k > 2 such that
    n is the m-th prime number, where m = (n+1)/k is also prime.

    This is equivalent to finding a prime m such that (p_m + 1) is divisible by m,
    and the quotient k = (p_m + 1) / m is greater than 2.
    """
    # We need to generate primes up to p_m. p_m is roughly m*ln(m).
    # To test m up to 30,000, p_m is roughly 30000*ln(30000) = 309,000.
    # A limit of 400,000 for our sieve is safe.
    LIMIT = 400000

    # Step 1: Generate primes using a Sieve of Eratosthenes.
    try:
        is_prime = [True] * (LIMIT + 1)
        is_prime[0] = is_prime[1] = False
        for p in range(2, int(math.sqrt(LIMIT)) + 1):
            if is_prime[p]:
                for multiple in range(p * p, LIMIT + 1, p):
                    is_prime[multiple] = False
        
        primes = [p for p, is_p in enumerate(is_prime) if is_p]
        num_primes = len(primes)
    except MemoryError:
        print("Error: The memory limit is too low to generate the primes. Please try on a machine with more RAM.")
        return

    # Step 2: Iterate through primes m to find a solution.
    # m is the index of the prime n, and m must also be prime.
    for m in primes:
        # We checked manually that for m=2 and m=3, k=2. The problem requires k>2.
        if m <= 3:
            continue

        # To find p_m (the m-th prime), we need primes[m-1].
        # This requires our list of primes to have at least m elements.
        if m > num_primes:
            # We cannot find the m-th prime if m is larger than the number of primes found.
            break

        # Step 3: Get n = p_m. Our list is 0-indexed, so p_m is at index m-1.
        n = primes[m - 1]

        # Step 4: Check if (n+1) is divisible by m.
        if (n + 1) % m == 0:
            # Step 5: Calculate k and check if it's greater than 2.
            k = (n + 1) // m
            if k > 2:
                # Solution found. Print the results and the equation.
                print(f"Solution found for k = {k} and n = {n}.")
                print(f"Here, n ({n}) is the {m}-th prime number, and m ({m}) is also prime.")
                print(f"The final equation is ({n} + 1) / {k} = {m}.")
                return

    # Step 7: If the loop finishes, no solution was found in the given range.
    print("Does not exist")

if __name__ == '__main__':
    find_special_prime()