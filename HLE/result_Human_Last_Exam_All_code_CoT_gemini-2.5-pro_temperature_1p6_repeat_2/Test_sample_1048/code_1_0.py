def find_best_permutation(a, x):
    """
    Finds the permutation of a list of numbers 'a' that maximizes the final
    value of 'x' after a series of modulo operations.

    This problem is solved using dynamic programming, similar to the
    Traveling Salesperson Problem.
    """
    n = len(a)
    initial_x = x

    # dp[mask][i] = max value using numbers in 'mask' ending with a[i]
    dp = [[-1] * n for _ in range(1 << n)]
    
    # parent[mask][i] = j, where a[j] was the element used before a[i]
    # in the optimal sequence for the subset 'mask'.
    parent = [[-1] * n for _ in range(1 << n)]

    # Base cases: sequences of length 1
    for i in range(n):
        dp[1 << i][i] = initial_x % a[i]

    # Fill the DP table for masks of size 2 to n
    for mask in range(1, 1 << n):
        # We process masks in increasing order of size implicitly
        # since prev_mask is always smaller than mask.
        for i in range(n):
            # If a[i] is in the current subset (mask)
            if (mask >> i) & 1:
                # We consider a[i] as the last element.
                # The previous mask did not contain a[i].
                prev_mask = mask ^ (1 << i)
                if prev_mask == 0:
                    continue
                
                max_val = -1
                best_j = -1
                # Iterate over possible previous elements a[j]
                for j in range(n):
                    if (prev_mask >> j) & 1:
                        # Value from the previous state followed by % a[i]
                        val = dp[prev_mask][j] % a[i]
                        if val > max_val:
                            max_val = val
                            best_j = j
                
                dp[mask][i] = max_val
                parent[mask][i] = best_j

    # Find the best result from all possible permutations
    final_mask = (1 << n) - 1
    best_x = -1
    last_idx = -1
    for i in range(n):
        if dp[final_mask][i] > best_x:
            best_x = dp[final_mask][i]
            last_idx = i

    # Reconstruct the optimal permutation by backtracking
    p = []
    curr_mask = final_mask
    curr_idx = last_idx
    while curr_idx != -1:
        p.insert(0, a[curr_idx])
        prev_idx = parent[curr_mask][curr_idx]
        curr_mask ^= (1 << curr_idx)
        curr_idx = prev_idx

    # Print the results as a full equation
    print("The optimal sequence of operations is:")
    
    # Use a temporary variable to show the calculation step-by-step
    temp_x = initial_x
    equation_str = str(initial_x)
    for val in p:
        equation_str += f" mod {val}"
        temp_x %= val
    
    equation_str += f" = {temp_x}"
    print(equation_str)

# --- Example Usage ---
# The counter-example for the greedy algorithm.
# The optimal result is 20, not 10.
list_a = [60, 40, 30]
start_x = 100

find_best_permutation(list_a, start_x)
