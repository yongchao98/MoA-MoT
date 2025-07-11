import sys

# Set a higher recursion limit for the depth of the calculation.
sys.setrecursionlimit(2000)

def solve():
    """
    Solves the problem of finding the largest union of 20 antichains in [1, N].
    """
    N = 823564528378596
    K = 20
    PRIMES = [2, 3, 5]
    MEMO = {}

    def count_ge(limit, k, p_idx):
        """
        Counts the number of integers n <= limit such that:
        1. n has at least k prime factors (counted with multiplicity).
        2. All prime factors of n are >= PRIMES[p_idx].
        """
        state = (limit, k, p_idx)
        if state in MEMO:
            return MEMO[state]

        # Base case: If k is non-positive, any integer qualifies.
        # The number of positive integers <= limit is floor(limit).
        if k <= 0:
            return limit

        # Base case: If we've exhausted our list of primes {2, 3, 5}, we stop.
        # Any number with >= 20 prime factors, all >= 7, is larger than N.
        if p_idx >= len(PRIMES):
            return 0

        p = PRIMES[p_idx]

        # Optimization: If the smallest possible number we can form, p^k,
        # is already greater than the limit, no such number exists.
        # We check p > 1 to avoid potential log(1) errors in other contexts,
        # and handle large exponents safely.
        if p > 1 and k * p.bit_length() > limit.bit_length() + 2: # Fast check for p**k > limit
             if p**k > limit:
                return 0

        # The count is the sum of two cases:
        # 1. Numbers not divisible by p (their smallest prime factor is > p).
        #    We move to the next prime in our list.
        res = count_ge(limit, k, p_idx + 1)
        
        # 2. Numbers divisible by p.
        #    These are p * m, where m <= limit/p and has >= k-1 prime factors >= p.
        res += count_ge(limit // p, k - 1, p_idx)

        MEMO[state] = res
        return res

    # According to a conjecture by Erdös, the largest union of 20 antichains is the set
    # of integers n <= N with Omega(n) <= 19.
    # We calculate this by finding the total number of integers (N) and subtracting
    # the count of integers with Omega(n) >= 20.
    
    # Calculate the number of integers with at least 20 prime factors.
    num_with_at_least_20_factors = count_ge(N, K, 0)
    
    # The final result is the total number of integers minus the "bad" ones.
    result = N - num_with_at_least_20_factors
    
    # Print the final equation as requested.
    print(f"{N} - {num_with_at_least_20_factors} = {result}")

solve()
<<<823564527457897>>>