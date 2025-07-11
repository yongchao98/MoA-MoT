import sys

# Set a higher recursion limit for deep recursive calls
sys.setrecursionlimit(2000)

# Problem parameters
N = 823564528378596
K = 19
PRIMES = [2, 3, 5]

# Memoization table
memo = {}

def psi3(x):
    """
    Counts numbers <= x that are not divisible by 2, 3, or 5,
    using the principle of inclusion-exclusion.
    This corresponds to counting numbers whose least prime factor is >= 7 (or the number 1).
    """
    if x < 1:
        return 0
    # Using integer division
    return (x - x//2 - x//3 - x//5 + 
            x//6 + x//10 + x//15 - 
            x//30)

def count_numbers(limit, k, p_idx):
    """
    Recursively counts numbers n <= limit such that Omega(n) <= k
    and the least prime factor of n is >= PRIMES[p_idx].
    """
    if k < 0:
        return 0
    if limit < 1:
        return 0
    
    # Memoization key
    state = (limit, k, p_idx)
    if state in memo:
        return memo[state]

    # Base case for recursion: all small primes (2, 3, 5) have been handled.
    if p_idx >= len(PRIMES):
        # For p>=7, p^(k+1) > N >= limit for any k<=19.
        # This means the condition Omega(n) <= k is always met for numbers
        # n <= limit with least prime factor >= 7.
        # So we just count such numbers.
        return psi3(limit)
    
    p = PRIMES[p_idx]

    # Recursive step
    # Count numbers with least prime factor >= PRIMES[p_idx + 1]
    res = count_numbers(limit, k, p_idx + 1)
    
    # Add count of numbers with prime factor p
    # These are p*m, where lpf(m) >= p and Omega(m) <= k-1
    res += count_numbers(limit // p, k - 1, p_idx)

    memo[state] = res
    return res

if __name__ == '__main__':
    # The final answer is the total count of numbers x <= N with Omega(x) <= 19.
    # We start the recursion with the first prime (index 0).
    total_size = count_numbers(N, K, 0)
    
    # Final result as requested by the user prompt's format.
    print(total_size)
    print(f"The largest union of 20 antichains is composed of all integers x from 1 to {N} where the number of prime factors of x (counted with multiplicity), Ω(x), is less than or equal to 19.")
    print(f"The final calculated size of this set is: {total_size}")