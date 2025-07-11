def solve_modulo_permutation():
    """
    Given a list of n numbers `a` and another number `x`, this function finds a
    permutation `p` of the list such that, after sequentially applying the
    operation x = x mod a_{p_i}, the resulting `x` has the smallest absolute
    difference from the original `x`.

    The function then prints the equation representing this optimal sequence of
    operations.
    """
    # Example values. You can change these to test with other inputs.
    a = [3, 8, 12]
    x_initial = 20
    n = len(a)

    # dp[mask] will store the set of reachable values using the subset of `a` represented by mask.
    dp = [set() for _ in range(1 << n)]
    
    # parent[mask][value] stores the previous state (mask, value) and the index of the number
    # used to transition to the current state. This is for backtracking.
    parent = [{} for _ in range(1 << n)]

    # Initial state: with no numbers used (mask=0), the only reachable value is x_initial.
    dp[0] = {x_initial}

    # Iterate through all possible subsets of `a` (represented by masks).
    for mask in range(1 << n):
        if not dp[mask]:
            continue
        
        # Try to add each unused number to the current subset.
        for i in range(n):
            if not (mask & (1 << i)):  # If i-th element is not in the current subset
                new_mask = mask | (1 << i)
                
                # For each value reachable with the current subset, calculate the new value.
                for val in dp[mask]:
                    new_val = val % a[i]
                    
                    # Add the new value to the set for the new subset.
                    # We store the parent only if this new_val is seen for the first time
                    # for this new_mask. This is sufficient to reconstruct one valid path.
                    if new_val not in dp[new_mask]:
                        dp[new_mask].add(new_val)
                        parent[new_mask][new_val] = (mask, val, i)

    final_mask = (1 << n) - 1
    
    # Handle the case of an empty list `a`.
    if n == 0:
        print(f"{x_initial} = {x_initial}")
        return
        
    # Find the best final value in the set of all possible outcomes.
    best_x_final = -1
    min_diff = float('inf')

    # Sort the list of final values to ensure the result is deterministic.
    for val in sorted(list(dp[final_mask])):
        diff = abs(val - x_initial)
        if diff < min_diff:
            min_diff = diff
            best_x_final = val

    # Reconstruct the permutation that leads to the best final value.
    perm_indices = []
    curr_mask = final_mask
    curr_val = best_x_final
    
    # Backtrack from the final state to the initial state.
    while curr_mask > 0:
        prev_mask, _, index = parent[curr_mask][curr_val]
        perm_indices.append(index)
        curr_mask = prev_mask
        # We need to find the value in the parent state that generated `curr_val`.
        # Since multiple paths can lead to the same (mask, val) state,
        # we need to find a valid predecessor. The `parent` table stores one.
        _, curr_val, _ = parent[curr_mask | (1 << index)][curr_val]

    perm_indices.reverse()

    # Build the string representing the sequence of operations.
    equation_str = str(x_initial)
    for index in perm_indices:
        equation_str = f"({equation_str} % {a[index]})"

    # Print the final result, including the equation.
    print(f"For a = {a} and x = {x_initial}:")
    print(f"The best final value is {best_x_final}, with a difference of {min_diff}.")
    print("The equation for the optimal permutation is:")
    print(f"{best_x_final} = {equation_str}")

solve_modulo_permutation()