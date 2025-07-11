import math

def solve():
    """
    Calculates the size of the largest union of 20 antichains in [1, N]
    in the divisor poset.
    This is equivalent to counting the number of integers n <= N with Omega(n) <= 19.
    We calculate this by finding the total N and subtracting numbers with Omega(n) >= 20.
    """
    N = 823564528378596
    
    # For k >= 20, any number n <= N with Omega(n)=k must have small prime factors,
    # since p^20 > N for p >= 7.
    # Primes to consider for the tail calculation (k>=20).
    # 7^19 < N, but 7^20 > N. So any number with Omega(n)=20 must have a prime factor < 7.
    primes = [2, 3, 5, 7]
    
    memo = {}

    def count(limit, k, p_idx):
        """
        Counts numbers n <= limit with Omega(n) = k,
        where each prime factor is >= primes[p_idx].
        """
        if k == 0:
            return 1
        if p_idx >= len(primes):
            return 0
        
        state = (limit, k, p_idx)
        if state in memo:
            return memo[state]

        p = primes[p_idx]
        
        # Pruning: if the smallest possible number p^k exceeds the limit,
        # then no such number exists.
        # Use math.log to avoid large integer exponentiation.
        if k * math.log(p) > math.log(limit) + 1e-9: # Add epsilon for float precision
            memo[state] = 0
            return 0

        # Sum counts for numbers whose smallest prime factor is p_j, where p_j >= p
        total = 0
        for i in range(p_idx, len(primes)):
            pi = primes[i]
            if pi > limit:
                break
            
            # Optimization: if pi^k > limit, subsequent primes will also result in numbers > limit
            if k * math.log(pi) > math.log(limit) + 1e-9:
                break

            total += count(limit // pi, k - 1, i)
        
        memo[state] = total
        return total

    # The maximum Omega(n) for n <= N is floor(log2(N)).
    k_max = int(math.log(N, 2))
    
    tail_sum = 0
    for k in range(20, k_max + 1):
        num_with_k_factors = count(N, k, 0)
        if num_with_k_factors == 0 and k > 25: # Optimization: if count is 0, likely 0 for larger k
            break
        tail_sum += num_with_k_factors
        
    result = N - tail_sum
    print(result)

solve()