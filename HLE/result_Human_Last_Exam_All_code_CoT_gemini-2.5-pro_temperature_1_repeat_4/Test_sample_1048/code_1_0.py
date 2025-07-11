import collections

def solve():
    """
    This function solves the problem of finding the permutation that maximizes the final x.
    It reads the input n, x, and the list a.
    Then, it implements the dynamic programming solution described above.
    """
    try:
        # Read n and x
        n, x = map(int, input().split())
        # Read the list a
        a = list(map(int, input().split()))
    except (IOError, ValueError):
        # Use a sample case if input reading fails
        print("Using sample case: n=3, x=29, a=[10, 15, 16]")
        n = 3
        x = 29
        a = [10, 15, 16]

    # dp[mask] will store the set of reachable values using the subset of 'a' represented by 'mask'.
    dp = collections.defaultdict(set)
    dp[0] = {x}

    # Iterate through all masks from 0 to 2^n - 1
    for mask in range(1 << n):
        # For each reachable value 'val' with the current mask
        for val in dp[mask]:
            # Try to apply the next number 'a[i]'
            for i in range(n):
                # Check if the i-th number has not been used yet
                if not (mask & (1 << i)):
                    # The next mask includes the i-th number
                    next_mask = mask | (1 << i)
                    # Calculate the new value and add it to the set for the next_mask
                    dp[next_mask].add(val % a[i])

    # The final set of reachable values is in dp[(1 << n) - 1]
    final_values = dp[(1 << n) - 1]

    if not final_values:
        # This case happens if n=0, the final x is the original x
        best_x = x
    else:
        # We want to maximize the final x to minimize the absolute difference
        best_x = max(final_values)
    
    original_x = x
    diff = abs(original_x - best_x)
    
    # We need to find the permutation that leads to best_x.
    # We can backtrack from the final state.
    
    # {value: mask}
    parent = {best_x: (1 << n) - 1}
    path_map = {best_x: []}
    
    q = collections.deque([(best_x, (1 << n) - 1)])
    
    final_permutation_found = False
    final_p = []

    # Backtracking search using BFS
    visited_states = set([(best_x, (1 << n) - 1)])

    while q:
        curr_val, curr_mask = q.popleft()

        if curr_mask == 0:
            if curr_val == x:
                final_p = path_map[curr_val]
                final_permutation_found = True
                break
            continue

        for i in range(n):
            if (curr_mask & (1 << i)):
                prev_mask = curr_mask ^ (1 << i)
                # Check all values in the previous state
                for prev_val in dp[prev_mask]:
                    if prev_val % a[i] == curr_val:
                        state = (prev_val, prev_mask)
                        if state not in visited_states:
                            path_map[prev_val] = [a[i]] + path_map[curr_val]
                            q.append(state)
                            visited_states.add(state)

    print(f"Original x: {original_x}")
    print(f"List a: {a}")
    print(f"The best possible final x is: {best_x}")
    print(f"The minimum absolute difference is: {diff}")
    
    # Print the equation step by step
    if final_p:
        print("One permutation that achieves this result is:", final_p)
        temp_x = original_x
        equation_str = f"{original_x}"
        for num in final_p:
            new_x = temp_x % num
            print(f"{temp_x} mod {num} = {new_x}")
            temp_x = new_x
    else:
        # This case for n=0
        print("No operations were performed.")


solve()