import math

def find_prime_solution():
    """
    Finds the smallest integer k > 2 and a prime number n such that n is the
    m-th prime number, where m = (n+1)/k is also prime.
    """
    primes = [2]
    num_to_check = 3
    # m_idx is the 0-based index for iterating through the 'primes' list.
    # The prime m will be primes[m_idx].
    m_idx = 0
    
    while True:
        # m is the (m_idx+1)-th prime number.
        m = primes[m_idx]

        # We need to find the m-th prime, n. To do this, our 'primes' list
        # must have at least m elements. The m-th prime is at index m-1.
        while len(primes) < m:
            is_prime = True
            # We check for primality by trial division against known primes up to sqrt(num).
            limit = int(math.sqrt(num_to_check))
            for p in primes:
                if p > limit:
                    break
                if num_to_check % p == 0:
                    is_prime = False
                    break
            if is_prime:
                primes.append(num_to_check)
            
            num_to_check += 2 # Move to the next odd number.
        
        # Now, the 'primes' list has at least 'm' prime numbers.
        # The m-th prime (n) is at index (m-1) in our 0-indexed list.
        n = primes[m - 1]

        # Now we check our main condition: k = (p_m + 1) / m must be an integer > 2.
        # This is equivalent to n = p_m, and k = (n + 1) / m.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            # Check if this k is the solution.
            if k > 2:
                print(f"Found the smallest values for k and n:")
                print(f"k = {k}")
                print(f"n = {n}")
                print("\nVerification:")
                print(f"The condition is that for k={k}, the number m = (n+1)/k must be prime, and n must be the m-th prime.")
                print(f"Let's calculate m: ({n} + 1) / {k} = {(n + 1) // k}")
                print(f"The value of m is {m}. Is it prime? Yes, it's the {m_idx + 1}-th prime.")
                print(f"Is n={n} the {m}-th prime? Yes, by our search method, it is.")
                break # Exit the loop as we've found the smallest solution.
        
        # Move to the next prime m for the next iteration.
        m_idx += 1

find_prime_solution()
<<<4, 43>>>