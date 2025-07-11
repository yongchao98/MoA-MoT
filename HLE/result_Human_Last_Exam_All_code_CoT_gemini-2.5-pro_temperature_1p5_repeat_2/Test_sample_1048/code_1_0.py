import math

def solve_modulo_permutation(a, x):
    """
    Finds the permutation of list 'a' that minimizes the absolute difference
    between the final and initial 'x' after sequential modulo operations.

    The method uses dynamic programming on subsets of 'a'.
    """
    if not a:
        print(f"The list is empty. Final x is {x}.")
        return

    n = len(a)
    x_initial = x

    # Find the minimum element 'm' and create A_prime (a without one m)
    min_val = min(a)
    min_idx = a.index(min_val)
    
    a_prime = a[:min_idx] + a[min_idx+1:]
    
    # dp[mask] = {value: (prev_mask, last_number_used)}
    # This stores backtracking info to reconstruct the permutation.
    dp = {0: {x_initial: None}}
    
    # Use a list of masks ordered by size to ensure correctness
    masks = sorted(range(1, 1 << len(a_prime)), key=bin)

    for mask in masks:
        dp[mask] = {}
        for i in range(len(a_prime)):
            if (mask >> i) & 1:
                prev_mask = mask ^ (1 << i)
                num = a_prime[i]
                if prev_mask in dp:
                    for val, _ in dp[prev_mask].items():
                        new_val = val % num
                        # Store path: new_val reached from val via num
                        dp[mask][new_val] = (prev_mask, val, num)

    # Find the best final x
    best_final_x = -1
    min_diff = float('inf')
    best_path_info = None

    all_reachable_inter_values = {}
    for mask, values in dp.items():
        for val in values:
            if val not in all_reachable_inter_values:
                 all_reachable_inter_values[val] = mask

    for inter_val, mask in all_reachable_inter_values.items():
        final_x = inter_val % min_val
        diff = abs(final_x - x_initial)

        if diff < min_diff:
            min_diff = diff
            best_final_x = final_x
            best_path_info = (mask, inter_val)
        elif diff == min_diff:
            if best_final_x == -1 or final_x < best_final_x:
                best_final_x = final_x
                best_path_info = (mask, inter_val)

    # Reconstruct the permutation prefix
    prefix = []
    if best_path_info:
        mask, val = best_path_info
        while mask > 0:
            prev_mask, prev_val, num = dp[mask][val]
            prefix.append(num)
            mask, val = prev_mask, prev_val
    prefix.reverse()

    # Construct the full permutation
    used_indices_in_aprime = set()
    temp_prefix = list(prefix)
    temp_aprime = list(a_prime)
    
    # Handle duplicates correctly when building the remainder
    # This part is to ensure we don't reuse elements for the remainder
    # that were already part of the optimal prefix
    final_prefix = []
    temp_a = list(a)
    temp_a.remove(min_val) # remove the m that will be used after prefix
    
    for p_val in prefix:
        final_prefix.append(p_val)
        temp_a.remove(p_val)
        
    permutation = final_prefix + [min_val] + temp_a

    # Print the equation
    equation = [str(x_initial)]
    current_val = x_initial
    for num in permutation:
        equation.append(f"% {num}")
        current_val %= num
    
    equation.append(f"= {current_val}")
    print(" ".join(equation))


# --- Example Usage ---
# You can change these values to test with other inputs.
# Note: This algorithm is exponential. It will be slow for n > 15.
a_list = [10, 4, 12]
x_val = 100
solve_modulo_permutation(a_list, x_val)

a_list_2 = [7, 4]
x_val_2 = 15
solve_modulo_permutation(a_list_2, x_val_2)
