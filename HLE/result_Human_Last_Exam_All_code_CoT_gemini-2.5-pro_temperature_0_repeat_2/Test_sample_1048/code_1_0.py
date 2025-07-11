import collections

def solve():
    """
    Solves the modulo permutation problem using dynamic programming on subsets.
    """
    try:
        line1 = input().split()
        if len(line1) < 2:
            print("Invalid input format. Please provide n, x, and the list a.")
            return
        n = int(line1[0])
        x = int(line1[1])
        a = list(map(int, input().split()))
        if len(a) != n:
            print(f"Input list size {len(a)} does not match n={n}.")
            return
    except (ValueError, IndexError):
        print("Invalid input. Please enter numbers only.")
        return

    # dp[mask] stores the set of achievable remainders for the subset `mask`.
    dp = collections.defaultdict(set)
    # parent[mask][remainder] stores how we achieved this state, for backtracking.
    # Format: (previous_mask, previous_value, number_used)
    parent = collections.defaultdict(dict)

    # Base case: with no numbers used, the only value is the original x.
    dp[0] = {x}

    # Iterate through all subsets of `a` from size 1 to n.
    for mask in range(1, 1 << n):
        # For each number `a[i]` in the current subset `mask`.
        for i in range(n):
            if (mask >> i) & 1:
                # This number `a[i]` is in the current set.
                # Consider it as the last number applied in a permutation of this subset.
                prev_mask = mask ^ (1 << i)
                
                # Calculate new remainders by applying `mod a[i]` to previous results.
                for val in dp[prev_mask]:
                    new_val = val % a[i]
                    dp[mask].add(new_val)
                    # Store parent pointer for backtracking.
                    # If multiple paths lead to the same state, any is fine.
                    parent[mask][new_val] = (prev_mask, val, a[i])

    # The set of all possible final outcomes.
    final_remainders = dp[(1 << n) - 1]

    if not final_remainders:
        best_x = x
    else:
        # Find the remainder that maximizes the final value (minimizes difference).
        best_x = max(final_remainders)

    # --- Backtracking to find the permutation ---
    perm = []
    curr_mask = (1 << n) - 1
    curr_val = best_x
    
    # Reconstruct the permutation from the parent pointers.
    while curr_mask > 0:
        prev_mask, prev_val, num_used = parent[curr_mask][curr_val]
        perm.append(num_used)
        curr_mask = prev_mask
        curr_val = prev_val
        
    perm.reverse()
    
    # --- Printing the output ---
    print(f"To minimize the absolute difference, we want to maximize the final remainder.")
    print(f"The best possible final remainder is: {best_x}")
    if perm:
        print(f"An optimal permutation is: {perm}")
    print("\nThe sequence of operations is:")
    
    current_x = x
    print(f"x = {current_x}")
    for num in perm:
        next_x = current_x % num
        print(f"x = {current_x} mod {num} = {next_x}")
        current_x = next_x
    
    print(f"\nFinal x = {current_x}")
    print(f"Original x = {x}")
    print(f"Absolute difference = |{current_x} - {x}| = {abs(current_x - x)}")

solve()