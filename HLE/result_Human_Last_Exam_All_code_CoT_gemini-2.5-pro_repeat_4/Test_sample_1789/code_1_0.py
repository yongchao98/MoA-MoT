import sys

# It's recommended to run this script with PyPy for better performance.
# The recursion limit is increased to handle deep recursive calls.
sys.setrecursionlimit(30000)

N = 823564528378596
K = 19

# The core of the solution is a recursive function (with memoization) to count
# the numbers with a certain number of prime factors.
# Let dp(x, k, i) be the number of integers n <= x such that Omega(n) <= k
# and the smallest prime factor of n is greater than or equal to the i-th prime, p_i.
# The recurrence relation is:
# dp(x, k, i) = dp(x, k, i+1)      (counts numbers whose smallest prime factor is > p_i)
#             + dp(x/p_i, k-1, i) (counts numbers for which p_i is a factor)
# The base cases for the recursion are:
# - If k < 0, the count is 0.
# - If x < 1, there are no numbers to count, so the result is 0.
# - If the current prime p_i > x, the only integer n <= x that can satisfy the
#   lpf(n) >= p_i condition is n=1. For n=1, Omega(1)=0, which is <= k. So the count is 1.

# A large number of primes may be needed for the recursion. We pre-compute them
# up to a reasonable limit. The algorithm's performance depends heavily on the
# number of states (x, k, i) that need to be computed and stored.
PRIME_LIMIT = 400000
primes = []
try:
    is_prime = [True] * (PRIME_LIMIT + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, PRIME_LIMIT + 1):
        if is_prime[p]:
            primes.append(p)
            for multiple in range(p * p, PRIME_LIMIT + 1, p):
                is_prime[multiple] = False
except MemoryError:
    print("MemoryError: Cannot allocate memory for the prime sieve. Try a smaller PRIME_LIMIT.")
    sys.exit(1)

memo = {}

def count_le_k_omega(x, k, i):
    """
    Counts the number of integers n <= x such that Omega(n) <= k
    and the smallest prime factor of n is >= primes[i].
    """
    if k < 0:
        return 0
    if x < 1:
        return 0
    
    # Base case for recursion on i: if the current prime is larger than x,
    # the only integer n <= x that can have a smallest prime factor >= primes[i] is n=1.
    # Omega(1) = 0, which is always <= k (since k>=0 here).
    if i >= len(primes) or primes[i] > x:
        return 1
        
    # Memoization to store and retrieve results for previously computed states.
    key = (x, k, i)
    if key in memo:
        return memo[key]
    
    # A crucial base case for the recursion on k.
    # If k=0, we are counting numbers with Omega(n)<=0, which is only n=1.
    # Since lpf(1) is effectively infinite, it's >= primes[i]. So, the count is 1.
    if k == 0:
        return 1

    p_i = primes[i]
    
    # The recurrence relation.
    res = count_le_k_omega(x, k, i + 1) + count_le_k_omega(x // p_i, k - 1, i)
    
    memo[key] = res
    return res

# The final result is obtained by calling the function with the initial parameters N, K,
# and starting with the first prime (index 0).
result = count_le_k_omega(N, K, 0)
print(result)