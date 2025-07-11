import sys

def solve_knapsack_multiple():
    """
    Solves the Multiple Knapsack Problem with unique item usage using dynamic programming.
    """
    # Knapsack Problem Details
    item_values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    item_weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    knapsack_capacities = [40, 120, 200]

    num_items = len(item_values)
    knapsack_capacities.sort()  # Sort for consistent indexing
    C0, C1, C2 = knapsack_capacities

    # Increase recursion limit for deep DP table access during backtracking, if needed.
    # This is a precaution for potentially deep recursion paths in some Python environments.
    # sys.setrecursionlimit(num_items + C0 + C1 + C2)
    
    # Initialize DP table: dp[i][c0][c1][c2]
    # Stores the max value for the first `i` items with capacities c0, c1, c2.
    # Using lists for simplicity, though numpy would be more memory-efficient.
    try:
        dp = [[[[0 for _ in range(C2 + 1)] for _ in range(C1 + 1)] for _ in range(C0 + 1)] for _ in range(num_items + 1)]
    except MemoryError:
        print("Error: The DP table is too large to fit in memory.")
        print("The required size is roughly (num_items * C0 * C1 * C2) integers.")
        return

    # Fill DP table
    for i in range(1, num_items + 1):
        item_index = i - 1
        v = item_values[item_index]
        w = item_weights[item_index]

        for c0 in range(C0 + 1):
            for c1 in range(C1 + 1):
                for c2 in range(C2 + 1):
                    # Option 1: Don't take the item
                    val_not_taken = dp[i-1][c0][c1][c2]

                    # Option 2: Place in knapsack 0
                    val_k0 = -1
                    if c0 >= w:
                        val_k0 = v + dp[i-1][c0 - w][c1][c2]

                    # Option 3: Place in knapsack 1
                    val_k1 = -1
                    if c1 >= w:
                        val_k1 = v + dp[i-1][c0][c1 - w][c2]

                    # Option 4: Place in knapsack 2
                    val_k2 = -1
                    if c2 >= w:
                        val_k2 = v + dp[i-1][c0][c1][c2 - w]
                    
                    # Store the maximum value among all options
                    dp[i][c0][c1][c2] = max(val_not_taken, val_k0, val_k1, val_k2)

    # Maximum value is in the final cell
    max_value = dp[num_items][C0][C1][C2]

    # Backtrack to find the items that form the solution
    selected_values = []
    c0, c1, c2 = C0, C1, C2
    for i in range(num_items, 0, -1):
        item_index = i - 1
        v = item_values[item_index]
        w = item_weights[item_index]
        
        # Check if the item was taken by comparing with the state without it
        if dp[i][c0][c1][c2] == dp[i-1][c0][c1][c2]:
            continue  # Item was not taken

        # If taken, determine which knapsack it was placed into
        if c0 >= w and dp[i][c0][c1][c2] == v + dp[i-1][c0 - w][c1][c2]:
            selected_values.append(v)
            c0 -= w
        elif c1 >= w and dp[i][c0][c1][c2] == v + dp[i-1][c0][c1 - w][c2]:
            selected_values.append(v)
            c1 -= w
        elif c2 >= w and dp[i][c0][c1][c2] == v + dp[i-1][c0][c1][c2 - w]:
            selected_values.append(v)
            c2 -= w

    # Print the final results
    print(f"Maximum possible total value: {max_value}")

    # Format and print the equation for the total value
    selected_values.sort(reverse=True)
    equation_str = " + ".join(map(str, selected_values))
    print(f"The values of the selected items are: {equation_str}")
    print(f"Final equation: {equation_str} = {max_value}")

if __name__ == '__main__':
    solve_knapsack_multiple()
<<<694>>>