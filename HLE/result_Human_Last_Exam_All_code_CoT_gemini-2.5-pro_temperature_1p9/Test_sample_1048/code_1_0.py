import math

def solve_modulo_permutation(a, x):
    """
    Finds the permutation of list 'a' that maximizes the final value of x
    after sequential a_i.

    The function uses dynamic programming to explore all possible outcomes.
    - dp[mask] stores a set of reachable values using the elements of 'a' indicated by 'mask'.
    - parent[(mask, val)] stores the previous state (prev_val, index_of_a) that led to 'val' with 'mask'.
    This allows reconstruction of the path that leads to the optimal final value.
    """
    n = len(a)
    initial_x = x

    # dp[mask] = set of reachable values
    dp = {0: {initial_x}}
    # parent[(mask, val)] = (previous_val, index_of_a_used)
    parent = {}

    for mask in range(1 << n):
        if mask not in dp:
            continue
        
        for val in dp[mask]:
            for i in range(n):
                # if i-th element is not used yet
                if not (mask & (1 << i)):
                    new_mask = mask | (1 << i)
                    new_val = val % a[i]
                    
                    if new_mask not in dp:
                        dp[new_mask] = set()

                    # Add the new state and store its parent for backtracking
                    if new_val not in dp[new_mask]:
                        dp[new_mask].add(new_val)
                        parent[(new_mask, new_val)] = (val, i)

    final_mask = (1 << n) - 1
    if final_mask not in dp or not dp[final_mask]:
        print("No solution found.")
        return

    # Find the best final x
    best_final_x = -1
    for val in dp[final_mask]:
        if val > best_final_x:
            best_final_x = val

    # Backtrack to find the path
    path = []
    curr_val = best_final_x
    curr_mask = final_mask
    
    for _ in range(n):
        prev_val, mod_idx = parent[(curr_mask, curr_val)]
        
        # Prepend to path to get the correct order
        path.insert(0, {
            'prev_val': prev_val,
            'mod': a[mod_idx],
            'new_val': curr_val
        })
        
        curr_val = prev_val
        curr_mask &= ~(1 << mod_idx)
    
    # Print the equation steps
    print(f"Initial x = {initial_x}")
    print("Optimal sequence of operations:")
    current_x_in_eq = initial_x
    final_eq_str = str(initial_x)
    for step in path:
        print(f"{step['prev_val']} % {step['mod']} = {step['new_val']}")
        # Append to the string for the full one-line equation
        final_eq_str = f"({final_eq_str} % {step['mod']})"

    final_permutation = [step['mod'] for step in path]
    print(f"\nThe optimal permutation is: {final_permutation}")
    print(f"The best resulting x is: {best_final_x}")
    print(f"The smallest absolute difference is: {abs(initial_x - best_final_x)}")

# Example usage:
# List a, and number x
a = [12, 8, 5]
x = 20
solve_modulo_permutation(a, x)

print("\n" + "="*20 + "\n")

a = [10, 3, 15]
x = 100
solve_modulo_permutation(a, x)
