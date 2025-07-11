import sys

def solve():
    """
    Finds the permutation of a list 'a' that maximizes the final value of 'x'
    after a series of modulo operations, and prints the sequence of operations.
    """
    # Read inputs
    # Example inputs:
    # a_str = "13 8 7 25"
    # x_initial_str = "123"
    try:
        a_str = input("Enter the list of numbers (space-separated): ")
        x_initial_str = input("Enter the initial number x: ")
        a = list(map(int, a_str.split()))
        x_initial = int(x_initial_str)
    except (ValueError, IndexError):
        print("Invalid input. Please enter numbers only.", file=sys.stderr)
        return

    n = len(a)
    if n == 0:
        print(f"The list is empty. The final x is {x_initial}.")
        return

    # dp[mask] is a dictionary mapping {value: (previous_value, last_index_used)}
    # It stores how a certain value was reached for a given subset of 'a'.
    dp = [{} for _ in range(1 << n)]
    dp[0] = {x_initial: (None, None)}

    # Iterate through all subsets of 'a' (represented by masks)
    for mask in range(1 << n):
        if not dp[mask]:
            continue
        
        # For each number not in the current subset
        for i in range(n):
            if not ((mask >> i) & 1):
                new_mask = mask | (1 << i)
                # For each value reachable with the current subset
                for val, _ in dp[mask].items():
                    new_val = val % a[i]
                    # Store the path to the new value. If a path already exists,
                    # it doesn't matter which one we keep, as they are equivalent
                    # for future calculations.
                    if new_val not in dp[new_mask]:
                        dp[new_mask][new_val] = (val, i)

    final_mask = (1 << n) - 1
    if not dp[final_mask]:
        print(f"No permutation can be applied. Final x is {x_initial}.")
        # This case happens if the list 'a' is empty, handled above.
        return

    # Backtrack to find the permutation
    best_final_x = max(dp[final_mask].keys())
    
    permutation = []
    current_val = best_final_x
    current_mask = final_mask

    while current_mask > 0:
        prev_val, prev_idx = dp[current_mask][current_val]
        permutation.append(a[prev_idx])
        current_val = prev_val
        current_mask ^= (1 << prev_idx)

    permutation.reverse()

    # Print the equation step by step
    print("\nOptimal sequence of operations:")
    x = x_initial
    print(f"Initial x = {x}")
    for num in permutation:
        next_x = x % num
        print(f"x = {x} % {num} = {next_x}")
        x = next_x
    
    print(f"\nThe original x was {x_initial}.")
    print(f"The best resulting x is {best_final_x}.")
    print(f"The absolute difference is |{best_final_x} - {x_initial}| = {abs(best_final_x - x_initial)}.")


solve()