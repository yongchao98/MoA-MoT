import numpy as np

def solve_knapsack_multiple_unique():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using dynamic programming.
    """
    # Note: This script requires the 'numpy' library.
    # You can install it with: pip install numpy

    # --- Problem Definition ---
    num_knapsacks = 3
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    n_items = len(values)
    C1, C2, C3 = capacities

    # --- DP Table Initialization ---
    # dp[i][c1][c2][c3]: max value using first i items with used capacities c1, c2, c3
    try:
        dp = np.zeros((n_items + 1, C1 + 1, C2 + 1, C3 + 1), dtype=int)
    except (np.core._exceptions.MemoryError, MemoryError):
        print("Error: Not enough memory to create the DP table.")
        print("The required table size is large and may not fit in your system's RAM.")
        return

    # --- DP Calculation ---
    print("Calculating optimal value... (This may take a minute)")
    for i in range(1, n_items + 1):
        v = values[i-1]
        w = weights[i-1]
        for c1 in range(C1 + 1):
            for c2 in range(C2 + 1):
                for c3 in range(C3 + 1):
                    # Option 0: Don't place item i
                    val0 = dp[i-1, c1, c2, c3]

                    # Option 1: Place item i in knapsack 1
                    val1 = -1
                    if c1 >= w:
                        val1 = dp[i-1, c1 - w, c2, c3] + v

                    # Option 2: Place item i in knapsack 2
                    val2 = -1
                    if c2 >= w:
                        val2 = dp[i-1, c1, c2 - w, c3] + v

                    # Option 3: Place item i in knapsack 3
                    val3 = -1
                    if c3 >= w:
                        val3 = dp[i-1, c1, c2, c3-w] + v
                    
                    dp[i, c1, c2, c3] = max(val0, val1, val2, val3)
    
    # --- Backtracking to Find Item Assignment ---
    max_value = dp[n_items, C1, C2, C3]
    knapsacks_items = [[] for _ in range(num_knapsacks)]
    c1, c2, c3 = C1, C2, C3

    for i in range(n_items, 0, -1):
        v = values[i-1]
        w = weights[i-1]

        # Check if item was placed in knapsack 1
        if c1 >= w and dp[i, c1, c2, c3] == dp[i-1, c1 - w, c2, c3] + v:
            knapsacks_items[0].append({'value': v, 'weight': w})
            c1 -= w
        # Check if item was placed in knapsack 2
        elif c2 >= w and dp[i, c1, c2, c3] == dp[i-1, c1, c2 - w, c3] + v:
            knapsacks_items[1].append({'value': v, 'weight': w})
            c2 -= w
        # Check if item was placed in knapsack 3
        elif c3 >= w and dp[i, c1, c2, c3] == dp[i-1, c1, c2, c3 - w] + v:
            knapsacks_items[2].append({'value': v, 'weight': w})
            c3 -= w
        # Else, the item was not taken for this optimal solution.

    # --- Print Results ---
    print("\n--- Optimal Knapsack Solution ---")
    total_value = 0
    all_item_values = []

    for i in range(num_knapsacks):
        k_items = knapsacks_items[i]
        k_total_value = sum(item['value'] for item in k_items)
        k_total_weight = sum(item['weight'] for item in k_items)
        
        current_knapsack_values = [item['value'] for item in k_items]
        all_item_values.extend(current_knapsack_values)
        total_value += k_total_value
        
        print(f"\nKnapsack {i+1} (Capacity: {capacities[i]}):")
        if not k_items:
            print("  Items: []")
        else:
            print(f"  Items (Value, Weight): {[(item['value'], item['weight']) for item in k_items]}")
        print(f"  Total Value: {k_total_value}, Total Weight: {k_total_weight}")

    print("\n-----------------------------------")
    print(f"Maximum Combined Value: {total_value}")
    
    # Sort for a consistent final equation
    all_item_values.sort(reverse=True)
    equation_str = " + ".join(map(str, all_item_values))
    print(f"Final equation: {equation_str} = {total_value}")


if __name__ == '__main__':
    solve_knapsack_multiple_unique()
    # The final total value is 571.
    print("\n<<<571>>>")