def solve():
    """
    Finds a permutation of a list of numbers that maximizes the final result
    of a series of modulo operations, and prints the full calculation.
    """
    # Inputs for the problem.
    # You can change these values to test with other cases.
    x_initial = 100
    a = [7, 13, 101, 20, 30]

    # Partition the list 'a' into numbers smaller/equal to x (A_small)
    # and numbers larger than x (A_large).
    A_small = [val for val in a if val <= x_initial]
    A_large = [val for val in a if val > x_initial]

    # If A_small is empty, no number can reduce x, so the result is x_initial.
    if not A_small:
        result_x = x_initial
        print(f"{x_initial}", end="")
        # The order of A_large elements doesn't matter, sort for deterministic output.
        for val in sorted(A_large):
            print(f" % {val}", end="")
        print(f" = {result_x}")
        return

    m = len(A_small)
    # dp[mask] will store the max result for the subset represented by the mask.
    dp = [-1] * (1 << m)
    # parent[mask] stores the index of the last element used to achieve dp[mask].
    parent = [-1] * (1 << m)
    
    # Base case: before any operations, the value is x_initial.
    dp[0] = x_initial

    # Iterate through all subsets of A_small, represented by masks.
    for mask in range(1, 1 << m):
        # For each element in the current subset...
        for i in range(m):
            if (mask >> i) & 1:  # Check if i-th element is in the mask
                # The mask for the subproblem (without the i-th element).
                prev_mask = mask ^ (1 << i)
                
                # Calculate the result if A_small[i] is the last element applied.
                val = dp[prev_mask] % A_small[i]
                
                # If this path gives a better (larger) result, update our tables.
                if val > dp[mask]:
                    dp[mask] = val
                    parent[mask] = i

    # Reconstruct the optimal permutation of A_small by backtracking from the final state.
    perm_small = []
    mask = (1 << m) - 1
    while mask > 0:
        last_element_index = parent[mask]
        perm_small.append(A_small[last_element_index])
        mask ^= (1 << last_element_index)
    perm_small.reverse()
    
    # The full permutation is A_large elements followed by the optimal ordering of A_small.
    full_permutation = sorted(A_large) + perm_small
    
    # Print the full equation step-by-step.
    current_x = x_initial
    print(f"{current_x}", end="")
    
    for val in full_permutation:
        print(f" % {val}", end="")
        current_x %= val
    
    print(f" = {current_x}")

solve()