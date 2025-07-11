import numpy as np

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using dynamic programming and backtracks to find the solution items.
    """
    # Knapsack problem inputs
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    n_items = len(values)
    n_knapsacks = len(capacities)

    # Initialize DP table. Using np.int32 for memory efficiency.
    # dp[i][w1][w2][w3] = max value for first i items with weights w1, w2, w3
    dp_shape = (n_items + 1, capacities[0] + 1, capacities[1] + 1, capacities[2] + 1)
    try:
        dp = np.zeros(dp_shape, dtype=np.int32)
    except MemoryError:
        print("Error: The DP table is too large for the available memory.")
        return

    # Fill the DP table
    for i in range(1, n_items + 1):
        v = values[i-1]
        w = weights[i-1]

        for w1 in range(capacities[0] + 1):
            for w2 in range(capacities[1] + 1):
                for w3 in range(capacities[2] + 1):
                    # Option 1: Do not include item 'i'
                    val_without_item = dp[i-1, w1, w2, w3]
                    
                    # Options to include item 'i' in one of the knapsacks
                    val_k1 = dp[i-1, w1-w, w2, w3] + v if w1 >= w else -1
                    val_k2 = dp[i-1, w1, w2-w, w3] + v if w2 >= w else -1
                    val_k3 = dp[i-1, w1, w2, w3-w] + v if w3 >= w else -1
                    
                    # Update DP table with the maximum value from all possibilities
                    dp[i, w1, w2, w3] = max(val_without_item, val_k1, val_k2, val_k3)

    # Find the maximum value from the last slice of the DP table
    final_dp_slice = dp[n_items]
    max_value = np.max(final_dp_slice)

    # Backtracking to find which items were chosen
    max_indices = np.where(final_dp_slice == max_value)
    w1, w2, w3 = max_indices[0][0], max_indices[1][0], max_indices[2][0]
    
    chosen_values = []
    
    # Iterate backwards from the last item
    for i in range(n_items, 0, -1):
        v = values[i-1]
        w = weights[i-1]
        current_val = dp[i, w1, w2, w3]

        # Check if the current item was placed in one of the knapsacks
        if w1 >= w and current_val == dp[i-1, w1-w, w2, w3] + v:
            chosen_values.append(v)
            w1 -= w  # Move to the state before this item was added
        elif w2 >= w and current_val == dp[i-1, w1, w2-w, w3] + v:
            chosen_values.append(v)
            w2 -= w
        elif w3 >= w and current_val == dp[i-1, w1, w2, w3-w] + v:
            chosen_values.append(v)
            w3 -= w
        # If none of the above, the item was not taken. The weight state for dp[i-1] is the same.

    chosen_values.reverse()

    # Print the final result
    print("Equation of chosen item values:")
    if chosen_values:
        equation_str = " + ".join(map(str, chosen_values))
        print(f"{equation_str} = {max_value}")
    else:
        print("No items selected = 0")
    
    print(f"\nMaximum possible total value:")
    print(max_value)


if __name__ == '__main__':
    solve_multiple_knapsack()
<<<737>>>