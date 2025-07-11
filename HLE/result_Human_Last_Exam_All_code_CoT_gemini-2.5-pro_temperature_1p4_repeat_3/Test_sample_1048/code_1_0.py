import collections

def solve():
    """
    Finds the optimal permutation of a list `a` to maximize the final value
    of `x` after a series of modulo operations.
    """
    # You can change the input list 'a' and the number 'x' here.
    a = [40, 30]
    x = 100

    n = len(a)
    initial_x = x

    # dp[mask] = set of achievable values using elements corresponding to the mask
    dp = collections.defaultdict(set)
    dp[0] = {x}
    
    # backtrack[mask][value] = (prev_mask, prev_value, number_used)
    # This dictionary stores predecessor information to reconstruct the path.
    backtrack = collections.defaultdict(dict)

    # Iterate through all possible subsets of `a`, represented by masks
    for mask in range(1 << n):
        if mask not in dp:
            continue
        # For each element not in the current subset
        for i in range(n):
            if not (mask & (1 << i)):
                new_mask = mask | (1 << i)
                # For each achievable value with the current subset
                for val in dp[mask]:
                    new_val = val % a[i]
                    # If we found a new achievable value for the next subset
                    if new_val not in dp[new_mask]:
                        dp[new_mask].add(new_val)
                        backtrack[new_mask][new_val] = (mask, val, a[i])

    final_mask = (1 << n) - 1
    best_x = -1

    if n == 0:
        best_x = x
    elif dp[final_mask]:
        # To minimize |initial_x - final_x|, we need to maximize final_x.
        best_x = max(dp[final_mask])

    # --- Outputting the results ---
    
    print(f"Given list a = {a} and x = {initial_x}")

    if best_x != -1 and n > 0:
        print("\nFinding the sequence of operations that yields the best result...")
        
        # Reconstruct the permutation and the calculation steps
        perm_steps = []
        curr_val = best_x
        curr_mask = final_mask
        
        while curr_mask > 0:
            prev_mask, prev_val, num_used = backtrack[curr_mask][curr_val]
            perm_steps.append((prev_val, num_used, curr_val))
            curr_mask = prev_mask
            curr_val = prev_val
        
        perm_steps.reverse()
        
        print("The final equation is built step-by-step as follows:")
        for prev_val, num_used, new_val in perm_steps:
            print(f"{prev_val} mod {num_used} = {new_val}")
        
        print(f"\nThe best resulting x is: {best_x}")
        print(f"The smallest absolute difference is: {abs(initial_x - best_x)}")

    elif n == 0:
        print(f"The list 'a' is empty. The final x is {x}.")
    else:
        print("No solution could be found.")

solve()