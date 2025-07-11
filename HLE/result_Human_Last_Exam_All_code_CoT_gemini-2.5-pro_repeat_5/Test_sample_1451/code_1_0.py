def solve():
    """
    Calculates the minimum number of trials T(n) for specified values of n.
    """
    ns = [2, 3, 1234, 6712]
    results = []

    # Find the maximum n to determine the range for our calculation
    max_n = max(ns)
    
    # memo_M will store M(k) = max(T_sub(j) for j from 0 to k)
    memo_M = {0: 0}

    # Iteratively compute M(k) for k from 1 to max_n
    for k in range(1, max_n + 1):
        # First, calculate T_sub(k)
        if k % 2 != 0:
            # If k is odd, the subproblem is solved with 0 additional weighings.
            t_sub_k = 0
        else:
            # If k is even, we perform k/2 weighings and face a new subproblem.
            # The worst case for that subproblem requires M(k/2) weighings.
            # M(k/2) would have been computed in a previous iteration.
            t_sub_k = k // 2 + memo_M[k // 2]
        
        # Next, calculate M(k)
        # M(k) is the maximum T_sub value seen so far up to k.
        # M(k-1) was computed in the previous iteration.
        memo_M[k] = max(memo_M[k - 1], t_sub_k)

    # Now that all M(k) values are computed, find T(n) for each required n
    for n in ns:
        # T(n) = n (initial weighings) + M(n) (worst-case subproblem weighings)
        final_T_n = n + memo_M[n]
        results.append(final_T_n)
        
    print(','.join(map(str, results)))

solve()