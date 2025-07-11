import math

def find_best_modulo_result(a, x_orig):
    """
    Finds a permutation of a list 'a' to apply sequentially to x = x % a[p_i]
    such that the final x is as close as possible to the original x.

    This function uses dynamic programming on subsets (bitmask DP).
    """
    n = len(a)

    if n == 0:
        print("The list 'a' is empty.")
        print(f"The final value is the original x: {x_orig}")
        print("Equation: x =", x_orig)
        return

    # dp[mask] is a dictionary mapping a reachable value to its path:
    # {value: (previous_value, index_of_last_operator)}
    dp = [{} for _ in range(1 << n)]
    dp[0] = {x_orig: None}  # Base case: mask=0, value=x_orig

    # Iterate through all subsets (masks)
    for mask in range(1 << n):
        # Skip unreachable masks
        if not dp[mask]:
            continue

        # For each value reachable with the current mask
        for val, _ in dp[mask].items():
            # Try to apply each unused element a[i]
            for i in range(n):
                # Check if the i-th element is NOT in the current mask
                if not (mask & (1 << i)):
                    new_mask = mask | (1 << i)
                    new_val = val % a[i]

                    # Store the path to this new value if it's the first time
                    # we are reaching this value with this mask.
                    if new_val not in dp[new_mask]:
                        dp[new_mask][new_val] = (val, i)

    # The final results are in the dp state for the all-ones mask
    final_values_map = dp[(1 << n) - 1]

    if not final_values_map:
        print("Error: Could not find a path that uses all numbers.")
        return

    # Find the best final x among all possible outcomes
    best_final_x = -1
    min_abs_diff = float('inf')

    for v in final_values_map.keys():
        diff = abs(v - x_orig)
        if diff < min_abs_diff:
            min_abs_diff = diff
            best_final_x = v
        elif diff == min_abs_diff:
            # If differences are equal, we prefer the larger final x,
            # as that corresponds to a smaller change from the original x.
            best_final_x = max(best_final_x, v)

    # Backtrack to reconstruct the optimal permutation
    permutation_indices = []
    current_val = best_final_x
    current_mask = (1 << n) - 1

    # We backtrack from the final state to the initial state
    for k in range(n, 0, -1):
        prev_val, last_idx = dp[current_mask][current_val]
        permutation_indices.append(last_idx)
        
        # Move to the previous state
        prev_mask = current_mask ^ (1 << last_idx)
        current_mask = prev_mask
        current_val = prev_val
    
    # The permutation was built backwards, so we reverse it
    permutation_indices.reverse()
    permutation_values = [a[i] for i in permutation_indices]

    # Print the results as requested
    print(f"Original x: {x_orig}")
    print(f"List a: {a}")
    print(f"The best final x is {best_final_x}, which has an absolute difference of {min_abs_diff} from the original.")
    print(f"Optimal permutation of 'a': {permutation_values}")
    
    print("\nCalculation steps:")
    x = x_orig
    equation_parts = [f"{x}"]
    for val in permutation_values:
        x_new = x % val
        equation_parts.append(f"mod {val} = {x_new}")
        x = x_new
    print(" -> ".join(equation_parts))


# Example usage:
# find_best_modulo_result(a=[6, 7], x_orig=100)
# find_best_modulo_result(a=[18, 12], x_orig=25)
# find_best_modulo_result(a=[10, 20, 50], x_orig=99)
