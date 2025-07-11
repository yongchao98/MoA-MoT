import math

def solve_modulo_permutation():
    """
    Finds the optimal permutation of a list `a` to maximize the final value of `x`
    after a series of sequential modulo operations.
    """
    # Example input. In a real application, this would be read from user input.
    a = [8, 9, 17]
    x_initial = 30
    n = len(a)

    print(f"Given numbers: {a}")
    print(f"Initial x: {x_initial}\n")

    # dp[mask] will store the set of reachable values using elements from `a`
    # corresponding to the set bits in the mask.
    dp = [set() for _ in range(1 << n)]
    dp[0] = {x_initial}

    # Fill the DP table by iterating through all masks.
    # The order of masks ensures that dp[prev_mask] is computed before dp[mask].
    for mask in range(1 << n):
        # Skip if the previous state is unreachable
        if not dp[mask]:
            continue
        
        for j in range(n):
            # If element a[j] is not in the current subset (mask), we can add it.
            if not (mask & (1 << j)):
                new_mask = mask | (1 << j)
                for val in dp[mask]:
                    dp[new_mask].add(val % a[j])

    # The set of all possible final outcomes is in the last dp entry.
    final_mask = (1 << n) - 1
    
    if not dp[final_mask]:
        best_x = x_initial
        p = []
    else:
        # Find the maximum possible final value of x.
        best_x = max(dp[final_mask])

        # Backtrack to find the permutation that leads to best_x.
        p = []
        current_val = best_x
        current_mask = final_mask

        for i in range(n, 0, -1):
            found_predecessor = False
            for j in range(n):
                # Check if j corresponds to an element in the current subset.
                if current_mask & (1 << j):
                    prev_mask = current_mask ^ (1 << j)
                    # Check if current_val could be obtained from a value in the predecessor set.
                    for prev_val in dp[prev_mask]:
                        if prev_val % a[j] == current_val:
                            # Found the predecessor, add a[j] to the permutation.
                            p.insert(0, a[j])
                            current_val = prev_val
                            current_mask = prev_mask
                            found_predecessor = True
                            break
                    if found_predecessor:
                        break
    
    # Print the results
    print(f"The permutation that yields the best result is: {p}")
    print(f"The best resulting x is: {best_x}")

    # Format and print the final equation.
    final_expr = str(x_initial)
    for val in p:
        final_expr = f"({final_expr} % {val})"
    
    print("\nThe calculation is:")
    print(f"{best_x} = {final_expr}")
    
    # Show the step-by-step calculation for clarity.
    print("\nStep-by-step:")
    current_x = x_initial
    print(f"Start with x = {current_x}")
    for i, val in enumerate(p):
        next_x = current_x % val
        print(f"Step {i+1}: x = {current_x} % {val} = {next_x}")
        current_x = next_x

solve_modulo_permutation()