import math

def solve_modulo_permutation(a, x):
    """
    Finds the permutation of list 'a' that minimizes the absolute difference
    between the final and original x after sequential modulo operations.

    Args:
        a (list[int]): A list of n numbers.
        x (int): The starting number.
    """
    n = len(a)
    if n == 0:
        print(f"The list 'a' is empty. The value of x remains {x}.")
        return

    original_x = x
    
    # dp[mask] will be a dictionary {value: (prev_value, last_index)}
    # It stores the reachable values for a given subset of 'a' (represented by mask)
    # and how to get there.
    dp = [{} for _ in range(1 << n)]
    
    # Base case: with no numbers used (mask=0), the only value is the original x.
    dp[0] = {original_x: (-1, -1)} # (-1, -1) is a sentinel for the start

    # Iterate through all subsets (masks)
    for mask in range(1 << n):
        if not dp[mask]:  # Skip masks that are not reachable
            continue
        
        # Try to apply the next unused number a[i]
        for i in range(n):
            # Check if the i-th number is not in the current subset
            if not ((mask >> i) & 1):
                next_mask = mask | (1 << i)
                # For each value reachable with the current subset
                for prev_val, _ in dp[mask].items():
                    new_val = prev_val % a[i]
                    # Store the new reachable value and the path to it
                    dp[next_mask][new_val] = (prev_val, i)

    # Find the best final value from all possible outcomes
    final_mask = (1 << n) - 1
    final_values = dp[final_mask]

    if not final_values:
        print("No solution found.")
        return

    best_val = -1
    min_diff = float('inf')

    for val in final_values.keys():
        diff = abs(val - original_x)
        if diff < min_diff:
            min_diff = diff
            best_val = val
        # If differences are equal, prefer the larger final value
        elif diff == min_diff and val > best_val:
            best_val = val

    # Reconstruct the permutation and the calculation steps
    path = []
    curr_val = best_val
    curr_mask = final_mask
    
    while curr_mask > 0:
        prev_val, i = dp[curr_mask][curr_val]
        path.append({'val': a[i], 'prev_x': prev_val, 'curr_x': curr_val})
        curr_mask ^= (1 << i)
        curr_val = prev_val
        
    path.reverse()

    # Print the step-by-step calculation
    print("Optimal sequence of operations:")
    for step in path:
        print(f"{step['prev_x']} mod {step['val']} = {step['curr_x']}")
    
    print(f"\nOriginal x: {original_x}")
    print(f"Final x: {best_val}")
    print(f"Absolute difference: {min_diff}")


# Example usage:
# Given values from the problem description
a = [6, 2, 7]
x = 20

# You can change the inputs here to test with other values
# a = [30, 40]
# x = 100

solve_modulo_permutation(a, x)