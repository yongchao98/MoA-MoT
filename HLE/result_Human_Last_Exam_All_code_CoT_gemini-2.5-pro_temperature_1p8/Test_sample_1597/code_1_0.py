def count_allowed_pairs():
    """
    This script calculates the number of ordered pairs of positive integers (a, b),
    where 1 <= a, b <= 1000, that are "allowed".
    
    An ordered pair (a, b) is allowed if for every primitive function f,
    f(ab) >= f(a)f(b). This condition is equivalent to:
    a = 1, or b = 1, or the sets of prime factors of a and b have a non-empty intersection.
    """
    
    N = 1000
    
    # A dictionary to store the computed prime factor sets to avoid re-computation.
    prime_factors_memo = {}

    # Step 1: Use a sieve to find the Smallest Prime Factor (spf) for numbers up to N.
    # This allows for efficient factorization.
    spf = list(range(N + 1))
    for i in range(2, int(N**0.5) + 1):
        if spf[i] == i:  # i is a prime number
            for j in range(i * i, N + 1, i):
                if spf[j] == j:  # if j's spf hasn't been set yet
                    spf[j] = i

    def get_prime_factors(n):
        """
        Returns the set of distinct prime factors of a number n
        using the pre-computed smallest prime factor array.
        """
        if n in prime_factors_memo:
            return prime_factors_memo[n]
        if n < 2:
            return set()
        
        factors = set()
        temp_n = n
        while temp_n > 1:
            p = spf[temp_n]
            factors.add(p)
            # Efficiently divide out all occurrences of this prime factor
            while temp_n % p == 0:
                temp_n //= p
        prime_factors_memo[n] = factors
        return factors

    # Pre-calculate all prime factor sets for numbers from 1 to N
    for i in range(1, N + 1):
        get_prime_factors(i)

    total_allowed_pairs = 0
    
    # Step 2: Iterate through 'a' from 1 to N and count allowed 'b's for each 'a'.
    for a in range(1, N + 1):
        # Case 1: a = 1
        # All pairs (1, b) are allowed.
        if a == 1:
            total_allowed_pairs += N
            continue
        
        # Case 2: a > 1
        # The allowed pairs are (a, 1) and pairs (a, b) where b > 1 and Omega(a) and Omega(b) have a common prime.
        # This is equivalent to 1 (for b=1) + number of b in [1, N] sharing a prime factor with a.
        
        p_factors_a = prime_factors_memo[a]
        p_factors = list(p_factors_a)
        k = len(p_factors)
        
        # Using Principle of Inclusion-Exclusion to count b's that share a prime factor with a.
        count_b_with_common_factor = 0
        
        # Iterate over all non-empty subsets of p_factors. (1 to 2^k - 1)
        for i in range(1, 1 << k):
            product_of_primes_in_subset = 1
            subset_size = 0
            for j in range(k):
                if (i >> j) & 1:
                    product_of_primes_in_subset *= p_factors[j]
                    subset_size += 1
            
            # The number of multiples of 'd' up to N is N // d.
            term = N // product_of_primes_in_subset
            
            # Add or subtract the term based on the subset size (Inclusion-Exclusion).
            if subset_size % 2 == 1: # Odd-sized subsets are added.
                count_b_with_common_factor += term
            else: # Even-sized subsets are subtracted.
                count_b_with_common_factor -= term
        
        # The number of allowed b's for this 'a' is 1 (for b=1) + the count calculated above.
        num_allowed_b = 1 + count_b_with_common_factor
        
        total_allowed_pairs += num_allowed_b
    
    # Print the final calculated number.
    print(f"The number of allowed ordered pairs (a,b) with 1 <= a,b <= 1000 is:")
    print(total_allowed_pairs)
    
    # Returning the value for the final answer format as well.
    return total_allowed_pairs

# Execute the function to get the answer.
final_answer = count_allowed_pairs()