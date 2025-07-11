import math

def solve():
    """
    Calculates the size of the largest union of 20 antichains in [1, N].
    This is equivalent to counting numbers n <= N with Omega(n) < 20.
    We compute the complement: count of numbers with Omega(n) >= 20 and subtract from N.
    """
    N = 823564528378596

    # For a number n with k prime factors, the smallest possible value is 2^k.
    # The largest prime factor p to consider must satisfy p^k <= N.
    # For the smallest k=20, p <= N^(1/20), which is less than 6.
    # The list of primes needed is very small.
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]

    memo = {}

    def count(k, limit, p_idx):
        """
        Counts numbers n <= limit with Omega(n) = k, where the smallest
        prime factor of n is >= primes[p_idx].
        """
        if k == 0:
            return 1
        
        state = (k, limit, p_idx)
        if state in memo:
            return memo[state]

        res = 0
        
        # Optimization: if the smallest possible number is already too large, prune the search.
        if math.log(primes[p_idx]) * k > math.log(limit):
             memo[state] = 0
             return 0

        for i in range(p_idx, len(primes)):
            p = primes[i]
            
            new_limit = limit // p
            if new_limit == 0:
                break
            
            res += count(k - 1, new_limit, i)
            
        memo[state] = res
        return res

    total_to_subtract = 0
    min_omega = 20
    max_omega = int(math.log2(N))  # This is 49 for the given N

    for k in range(min_omega, max_omega + 1):
        # Clear memo for each k to avoid large memory usage, as states for different k are independent.
        memo.clear()
        count_for_k = count(k, N, 0)
        total_to_subtract += count_for_k

    final_answer = N - total_to_subtract

    print(f"{N} - {total_to_subtract} = {final_answer}")
    print(f"\n<<< {final_answer} >>>")

solve()