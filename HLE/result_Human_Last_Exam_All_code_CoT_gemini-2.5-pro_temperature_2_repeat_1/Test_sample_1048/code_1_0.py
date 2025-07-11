import collections

def find_best_permutation(a, x_orig):
    """
    Finds the best permutation of `a` for sequential modulo operations on `x_orig`.

    This function uses dynamic programming over subsets of `a`.
    `dp[mask]` is a dictionary mapping achievable values to their backtracking information.
    dp[mask][value] = (prev_value, last_index)
    - `mask`: A bitmask representing the subset of `a` used.
    - `value`: An achievable remainder after the operations.
    - `prev_value`: The value before the last modulo operation.
    - `last_index`: The index in `a` of the number used in the last operation.
    """
    n = len(a)
    if n == 0:
        return [], x_orig

    # dp[mask] -> {value: (prev_value, last_element_index)}
    dp = collections.defaultdict(dict)
    dp[0] = {x_orig: (None, -1)}

    # Iterate through all subsets of a, represented by masks
    for mask in range(1, 1 << n):
        for i in range(n):
            if (mask >> i) & 1:  # If the i-th element is in the current subset
                prev_mask = mask ^ (1 << i)
                if prev_mask in dp:
                    for prev_val, _ in dp[prev_mask].items():
                        new_val = prev_val % a[i]
                        # Store how we got here. If new_val can be reached in multiple ways,
                        # we only need to store one path.
                        if new_val not in dp[mask]:
                            dp[mask][new_val] = (prev_val, i)

    final_mask = (1 << n) - 1
    # Check if the final state is reachable (it should be for n>0).
    if not dp[final_mask]:
      best_val = x_orig
    else:
        # Find the value in the final DP state that's closest to x_orig
        min_diff = float('inf')
        best_val = -1
        for val in dp[final_mask]:
            diff = abs(val - x_orig)
            if diff < min_diff:
                min_diff = diff
                best_val = val
            # Tie-breaking rule: prefer the larger remainder if differences are equal
            elif diff == min_diff and val > best_val:
                best_val = val

    # Backtrack to reconstruct the permutation
    p_indices = []
    current_mask = final_mask
    current_val = best_val
    while current_mask > 0:
        prev_val, last_idx = dp[current_mask][current_val]
        p_indices.append(last_idx)
        current_mask ^= (1 << last_idx)
        current_val = prev_val

    p_indices.reverse()
    permutation = [a[i] for i in p_indices]
    
    return permutation, best_val

# --- Example Execution ---

# You can change these values to test with other inputs
a_list = [100, 2, 3]
x_start = 12345

# Find the optimal permutation and the resulting x
permutation, final_x = find_best_permutation(a_list, x_start)

# Format the final equation for printing
equation_parts = [str(x_start)]
for num in permutation:
    equation_parts.append(f"% {num}")

equation_str = " ".join(equation_parts)
print(f"The best permutation results in the following equation:")
print(f"{equation_str} = {final_x}")
