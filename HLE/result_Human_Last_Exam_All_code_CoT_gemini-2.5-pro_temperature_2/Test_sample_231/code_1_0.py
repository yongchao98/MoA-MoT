import math

def find_k_and_n():
    """
    This function searches for a prime number n and an integer k > 2 such that:
    1. n is the p-th prime number, where p = (n+1)/k.
    2. p is also a prime number.

    This can be rewritten as: n = P_p and n+1 = k*p, where n and p are prime, and k > 2.
    
    The strategy is to iterate through candidate primes p, find the corresponding n=P_p,
    and then check if (n+1) is divisible by p to find an integer k > 2.
    
    Note: This search can be computationally intensive and may take a few minutes to run.
    """

    # A list to store generated primes. We start with the first prime.
    primes_list = [2]
    # We start checking for new primes from the next odd number.
    num_to_check = 3
    
    # We will check primes 'p' using their index in our list.
    p_candidate_idx = 0 

    # We set a search limit to avoid an infinite loop. 
    # This limit is on the index of the prime p, not p itself.
    # A limit of 3500 means we test the first 3500 primes as candidates for p.
    # The 3500th prime is 32443.
    search_idx_limit = 3500

    while p_candidate_idx < search_idx_limit:
        # Step 1: Get the current prime candidate 'p'.
        # We need to ensure primes_list is long enough to provide primes_list[p_candidate_idx].
        # This is implicitly handled as the list grows, but we need to be careful with the loop logic.
        
        # This while loop ensures our list of primes is large enough to get the next p candidate.
        while p_candidate_idx >= len(primes_list):
            is_prime_num = True
            sqrt_num = int(math.sqrt(num_to_check))
            for prime in primes_list:
                if prime > sqrt_num:
                    break
                if num_to_check % prime == 0:
                    is_prime_num = False
                    break
            if is_prime_num:
                primes_list.append(num_to_check)
            # Always check the next odd number.
            num_to_check += 2
            
        p = primes_list[p_candidate_idx]

        # Step 2: Find 'n', which is the p-th prime.
        # This means n is at index (p-1) in our 0-indexed list.
        # We need to extend our list of primes until it has at least p elements.
        while len(primes_list) < p:
            is_prime_num = True
            sqrt_num = int(math.sqrt(num_to_check))
            for prime in primes_list:
                if prime > sqrt_num:
                    break
                if num_to_check % prime == 0:
                    is_prime_num = False
                    break
            if is_prime_num:
                primes_list.append(num_to_check)
            num_to_check += 2

        n = primes_list[p - 1]

        # Step 3: Check if (n + 1) is divisible by p.
        if (n + 1) % p == 0:
            k = (n + 1) // p
            # Step 4: Check if k > 2.
            if k > 2:
                # We found a solution. Because we iterate on p in increasing order,
                # this finds the solution with the smallest n and k.
                print(f"{k}, {n}")
                return

        # Move to the next candidate for p.
        p_candidate_idx += 1

    # If the loop completes without finding a solution, we report that none was found.
    print("Does not exist")


if __name__ == '__main__':
    find_k_and_n()